library(DescTools)
# variables
max_exon_pos = -4
max_intron_pos = 8
window_length = 250
regex_donor <- paste0("(?=([A-Z]{", abs(max_exon_pos), "}G[TC]{1}[A-Z]{", (max_intron_pos-2), "}))")

# get ALT sequence from variant information
get_alt_seq <- function(region_start, var_start, ref, alt, ref_seq) {
  alt = gsub(pattern = paste0("(\\w{", var_start - region_start, "})\\w{", nchar(ref), "}(.*)"),
             replacement = paste0("\\1", alt, "\\2"),
             x = ref_seq, perl = T)
  return(alt)
}
  
#annotate whether variant is SNP, DEL or INS
variant_type <- function(ref, alt) {
  if (nchar(ref) == nchar(alt)) {vt = 'SNP'
  } else if (nchar(ref) > 1) {vt = 'DEL'
  } else if (nchar(alt) > 1) {vt = 'INS'
  }
  return(vt)
}

# find variant position start/end relative to authentic donor
find_dist_from_donor_start <- function(ADP, VS, VE, strand, variant_type) {
  dfds = NA
  sa <- ifelse(strand == '+', 1, -1)
  VS_sa <- ifelse(strand == '+', VS, VE)
  pa <- ifelse( (sa * VS_sa) >= (sa * ADP), 1, 0)
  
  dfds = sa * (VS_sa - ADP) + pa
  
  return(dfds)
}


# annotate CSS variant category
## category A: pre-existing cryptic, variant affects authentic 5'SS
## category B: variant-created cryptic, variant affects only cryptic 5'SS
## category C: variant-created cryptic overlapping with authentic 5'SS
## category D: pre-existing cryptic, variant does not affect authentic or cryptic 5'SS
'%!overlaps%' <- function(x,y)!('%overlaps%'(x,y))

find_cryptic_category <- function(ADP, VS, VE, DS, DE,CS, CE) {
  donor_range = c(DS, DE)
  variant_range = c(VS, VE)
  cryptic_range = c(CS, CE)
  
  if (variant_range %overlaps% donor_range & variant_range %!overlaps% cryptic_range) {
    category = 'A'
  } else if (variant_range %!overlaps% donor_range & variant_range %overlaps% cryptic_range) {
    category = 'B'
  } else if (variant_range %overlaps% donor_range & variant_range %overlaps% cryptic_range) {
    category = 'C'
  } else if (variant_range %!overlaps% donor_range & variant_range %!overlaps% cryptic_range) {
    category = 'D'
  }else (category = NA)
    
  return(category)
}


annotate_variant_category <- function(nested_cryptic_category) {
  cryp_cats <- sort(unlist(strsplit(nested_cryptic_category, split = '\\;')))
  if (identical(unique(cryp_cats),c('A'))) {
    var_cat = 'A'
  } else if (identical(unique(cryp_cats),c('B'))) {
    var_cat = 'B'
  } else if (identical(unique(cryp_cats),c('C'))) {
    var_cat = 'C'
  } else if (identical(unique(cryp_cats),c('A','C'))) {
    var_cat = 'C'
  } else (var_cat = NA)
  return(var_cat)
}

find_dist_from_donor_end <- function(ADP, VS, VE, strand, variant_type) {
  dfde = NA
  sa <- ifelse(strand == '+', 1, -1)
  VE_sa <- ifelse(strand == '+', VE, VS)
  pa <- ifelse( (sa * VE_sa) >= (sa * ADP), 1, 0)
  
  dfde = sa * (VE_sa - ADP) + pa
  
  return(dfde)
}



# find position in strings of ref and var authentic donor, accounting for strand and indels - to calculate NIF
find_ref_sp <- function(strand, ADP, RS, RE) {
  ref_sp <- NA
  strand_adj <- ifelse(strand == '+', 1, -1)
  RS_SA <- ifelse(strand == '+', RS, RE)
  ref_sp <- (strand_adj * (ADP - RS_SA))
  ref_sp <- ref_sp + 1
  return(ref_sp)
}

find_alt_sp <- function(strand, variant_type, ADP, VS, VE, RS, RE, ref, alt) {
  alt_sp <- NA
  strand_adj <- ifelse(strand == '+', 1, -1)
  RS_SA <- ifelse(strand == '+', RS, RE)
  VS_SA <- ifelse(strand == '+', VS, VE)
  VE_SA <- ifelse(strand == '+', VE, VS)
  indel_adj <- nchar(alt) - nchar(ref)
  
  alt_sp <- (strand_adj * (ADP - RS_SA))
  
  # adjustments for indels
  if ((variant_type == 'INS' & (strand_adj * VS_SA) < (strand_adj * ADP))|
      (variant_type == 'DEL' & (strand_adj *  VS_SA) < (strand_adj * ADP) & 
       (strand_adj *  VE_SA) < (strand_adj * ADP))) {
    alt_sp <- alt_sp + indel_adj
  } else if (variant_type == 'DEL' & (strand_adj *  VS_SA) < (strand_adj * ADP) & 
             (strand_adj *  VE_SA) >= (strand_adj * ADP) & strand == '+') {
    alt_sp <- (strand_adj * (VS_SA + 1 - RS_SA)) 
  } else if (variant_type == 'DEL' & (strand_adj *  VS_SA) < (strand_adj * ADP) & 
             (strand_adj *  VE_SA) >= (strand_adj * ADP) & strand == '-') {
    alt_sp <- (strand_adj * (VS_SA - RS_SA)) 
  }
  alt_sp <- alt_sp + 1
  return(alt_sp)
}




find_var_start_sp <- function(variant_type, strand, authentic_alt_sp, VS, VE, ADP,ref,alt) {
  var_start_sp = NA
  strand_adj <- ifelse(strand == '+', 1, -1)
  indel_adj <- nchar(alt) - nchar(ref)
  VS_SA <- ifelse(strand == '+', VS, VE)
  VE_SA <- ifelse(strand == '+', VE, VS)
  special_case <- ifelse(strand == '+', 1, 0)
  
  var_start_sp = authentic_alt_sp + strand_adj * (VS_SA - ADP)
  
  if (variant_type == 'INS' & (strand_adj * VS_SA) < (strand_adj * ADP)) {
    var_start_sp = var_start_sp - indel_adj
    
  } else if (variant_type == 'DEL' & (strand_adj * VS_SA) < (strand_adj * ADP) & 
             (strand_adj * VE_SA) < (strand_adj * ADP) ) {
    var_start_sp = authentic_alt_sp + strand_adj * (VE_SA - ADP)
    
  } else if (variant_type == 'DEL' & (strand_adj * VS_SA) < (strand_adj * ADP) &
             (strand_adj * VE_SA) >= (strand_adj * ADP)) {
    var_start_sp = authentic_alt_sp - special_case
    
  }
   return(var_start_sp) 
}

find_var_end_sp <- function(variant_type, strand, authentic_alt_sp, VS, VE, ADP,ref,alt) {
  var_end_sp = NA
  strand_adj <- ifelse(strand == '+', 1, -1)
  indel_adj <- nchar(alt) - nchar(ref)
  VS_SA <- ifelse(strand == '+', VS, VE)
  VE_SA <- ifelse(strand == '+', VE, VS)
  special_case <- ifelse(strand == '+', 1, 0)
  
  var_end_sp = authentic_alt_sp + strand_adj * (VS_SA - ADP)
  
  if (variant_type == 'INS' & (strand_adj * VS_SA) >= (strand_adj * ADP)) {
    var_end_sp = var_end_sp + indel_adj
    
  } else if (variant_type == 'DEL' & (strand_adj * VS_SA) < (strand_adj * ADP) & 
             (strand_adj * VE_SA) < (strand_adj * ADP) ) {
    var_end_sp = authentic_alt_sp + strand_adj * (VE_SA - ADP)
    
  } else if (variant_type == 'DEL' & (strand_adj * VS_SA) < (strand_adj * ADP) &
             (strand_adj * VE_SA) >= (strand_adj * ADP)) {
    var_end_sp = authentic_alt_sp - special_case
    
  }
  return(var_end_sp) 
}



# finding decoy donor genomic location
find_decoy_gloc <- function(variant_type, strand, ADP, decoy_distance, var_start, var_end, 
                            aa_sp, decoy_sp, var_start_sp, var_end_sp, ref, alt) {
  decoy_gloc = NA
  sa <- ifelse(strand == '+', 1, -1)
  indel_adj <- nchar(ref) - nchar(alt)
  del_sc_start <- ifelse(strand == '+', var_start + 1, var_end)
  del_sc_end <- ifelse(strand == '+', var_end + 1, var_start)
  
  decoy_gloc =  ADP + sa * decoy_distance
  
  if (variant_type == 'DEL' & (sa * var_start) < (sa *( ADP - 1)) & (sa * var_end) < (sa * ADP) & decoy_sp < var_start_sp) {
    decoy_gloc = decoy_gloc - sa * indel_adj
  } else if (variant_type == 'DEL' & (sa * var_start) > (sa * (ADP - 1)) & (sa * var_end) > (sa * ADP) & decoy_sp > var_start_sp) {
    decoy_gloc = decoy_gloc + sa * indel_adj
    
  } else if (variant_type == 'DEL' & var_start <= ADP - 1 & var_end >= ADP & decoy_sp < var_start_sp) {
    decoy_gloc = del_sc_start + sa * decoy_distance
  } else if (variant_type == 'DEL' & var_start <= ADP - 1 & var_end >= ADP & decoy_sp > var_start_sp) {
    decoy_gloc = del_sc_end + sa * decoy_distance
    
  } else if (variant_type == 'INS' & decoy_sp >= var_start_sp & decoy_sp <= var_end_sp) {
    decoy_gloc = var_start
  } else if (variant_type == 'INS' & (sa*var_start) < (sa*ADP) & decoy_sp < var_start_sp) {
    decoy_gloc = decoy_gloc - sa * indel_adj
  } else if (variant_type == 'INS' & (sa* var_start) >= (sa*ADP) & decoy_sp > var_end_sp) {
    decoy_gloc = decoy_gloc + sa * indel_adj
  }
  return(decoy_gloc)
}


var_dist_from_cryptic <- function(var_pos, cryptic_pos, strand) {
  dist = NA
  if (strand == '+') {
    
    dist <- var_pos - cryptic_pos
    
    if (var_pos >= cryptic_pos) {
      dist <- dist + 1 #intronic variants
    }
    
  } else if (strand == '-') {
    dist <- cryptic_pos - var_pos
    
    if (cryptic_pos >= var_pos) {
      dist <- dist + 1
    }
  }
  
  return(dist)
}

# find resultant intron and exon width after decoy / cryptic use
calculate_resultant_intron_width <- function(donor_string_position, alt_seq) {
  riw = NA 
  riw = nchar(alt_seq) - donor_string_position - 7 + 1
  return(riw)
}

calculate_resultant_exon_width <- function(donor_string_position) {
  rew = NA 
  rew = donor_string_position - 4 - 1
  return(rew)
}


# Functions for calculating decoy depletion
decoy_depletion_grouped <- function(simulation_df, obs_decoy_df, grouping_variable = NULL, decoy_dist_filter = 250) {
  if (is.null(grouping_variable)) {
    simulation_decoy_counts <- simulation_df[abs(decoy_loc) <= decoy_dist_filter][!is.na(decoy_loc), .(n = .N), 
                                                                                  by = c('decoy_loc', 
                                                                                         'simulation_no')][, .(simulated_decoys = mean(n)), 
                                                                                                           by = c('decoy_loc')][order(decoy_loc)]
    
    # summarise the number of decoys per position in hg19, compared to simulations
    obs_decoy_counts <- simulation_decoy_counts[obs_decoy_df[abs(decoy_loc) <= decoy_dist_filter][!is.na(decoy_loc), 
                                                                                                  .(observed_decoys= .N), 
                                                                                                  by = c('decoy_location','decoy_loc', 'frame')], 
                                                on = c('decoy_loc')][order(decoy_loc)]
    
    obs_decoy_counts[, obs_v_exp := observed_decoys / simulated_decoys]
    
    return(obs_decoy_counts)
  } else {
    # summarise the mean number of decoys per position across 15 simulations
    simulation_decoy_counts <- simulation_df[abs(decoy_loc) <= decoy_dist_filter][!is.na(decoy_loc), .(n = .N), 
                                                                                  by = c('decoy_loc', 
                                                                                         'simulation_no', 
                                                                                         grouping_variable)][, .(simulated_decoys = mean(n)), 
                                                                                                             by = c('decoy_loc', grouping_variable)][order(decoy_loc)]
    
    # summarise the number of decoys per position in hg19, compared to simulations
    obs_decoy_counts <- simulation_decoy_counts[obs_decoy_df[abs(decoy_loc) <= decoy_dist_filter][!is.na(decoy_loc), .(observed_decoys= .N), 
                                                                                                  by = c('decoy_location','decoy_loc',grouping_variable)], 
                                                on = c('decoy_loc', grouping_variable)][order(decoy_loc)]
    
    obs_decoy_counts[, obs_v_exp := observed_decoys / simulated_decoys]
    
    return(obs_decoy_counts)
  }
  
}


plot_decoy_depletion <- function(simulation_df, obs_decoy_df, grouping_variable = NULL, decoy_dist_filter = 250) {
  if (is.null(grouping_variable)) {
    plot <- decoy_depletion_grouped(simulation_df, obs_decoy_df, decoy_dist_filter = decoy_dist_filter) %>%
      ggplot(aes(x = decoy_loc, y = obs_v_exp)) +
      geom_point(size = 0.2) + facet_wrap(~decoy_location, scales = 'free_x', ncol = 3) + 
      geom_smooth(aes(color = decoy_location))
  } else {
    plot <- decoy_depletion_grouped(simulation_df, obs_decoy_df, grouping_variable, decoy_dist_filter) %>%
      ggplot(aes_string(x = 'decoy_loc', y = 'obs_v_exp', color = grouping_variable)) +
      geom_point(size = 0.2) + facet_wrap(~decoy_location, scales = 'free_x', ncol = 3) + 
      geom_smooth()
  }
  return(plot)
}

plot_decoy_depletion_detailed <- function(simulation_df, obs_decoy_df, grouping_variable = NULL, decoy_dist_filter = 250) {
  if (is.null(grouping_variable)) {
    plot <- pivot_longer(decoy_depletion_grouped(simulation_df, obs_decoy_df, decoy_dist_filter = decoy_dist_filter), 
                         cols = c('simulated_decoys','observed_decoys', 'obs_v_exp')) %>%
      mutate(name = factor(name, levels = c('simulated_decoys','observed_decoys', 'obs_v_exp'))) %>%
      ggplot(aes(x = decoy_loc, y = value)) +
      geom_point(size = 0.2) + facet_wrap(name~decoy_location, scales = 'free', ncol = 2)
  } else {
    plot <- pivot_longer(decoy_depletion_grouped(simulation_df, obs_decoy_df, grouping_variable, decoy_dist_filter), 
                         cols = c('simulated_decoys','observed_decoys', 'obs_v_exp')) %>%
      mutate(name = factor(name, levels = c('simulated_decoys','observed_decoys', 'obs_v_exp'))) %>%
      ggplot(aes_string(x = 'decoy_loc', y = 'value', color = grouping_variable)) +
      geom_point(size = 0.2) + facet_wrap(name~decoy_location, scales = 'free', ncol = 2)
  }
  return(plot)
}


