# All functions required for 1c_Extract_Metabolites_M0M1M2_using_Identification.R

# Function countElements that convert the chemical formula string to a recognized numerical format
countElements <- function(inputformula) {
  elements <- c("H", "C", "N", "O", "S", "P", "Si")
  counts <- rep(0, length(elements))
  
  hasSi <- grepl("Si", inputformula)
  
  for (i in 1:length(elements)) {
    element <- elements[i]
    
    if (element == "S") {
      modified_formula <- inputformula
      if (hasSi) {
        modified_formula <- gsub("Si\\d*", "", modified_formula)  # Remove "Si" and the value after it
      }
      
      pattern <- paste0(element, "\\d*")
      counts[i] <- sum(sapply(regmatches(modified_formula, gregexpr(pattern, modified_formula)), function(x) {
        if (length(x) > 0) {
          num <- gsub(element, "", x)
          if (num == "") num <- 1
          as.numeric(num)
        } else {
          0
        }
      }))
    } else {
      pattern <- paste0(element, "\\d*")
      counts[i] <- sum(sapply(regmatches(inputformula, gregexpr(pattern, inputformula)), function(x) {
        if (length(x) > 0) {
          num <- gsub(element, "", x)
          if (num == "") num <- 1
          as.numeric(num)
        } else {
          0
        }
      }))
    }
  }
  
  result <- data.frame(Element = elements, Count = counts)
  return(result)
}

# Function calculate_first_derivative and calculate_second_derivative to help determine peak cutoff, 
# the calculating method is mimicking the MS-DIAL method
calculate_first_derivative <- function(intCol) {
  total_slices <- length(intCol)
  # Initialize first derivative intensities vector
  fir_derivative <- numeric(total_slices)
  
  # Loop through each slice
  for (p in 1:total_slices) {
    if (p <= 2 || p >= total_slices - 1) {
      # For the first two and last two data points, set the first derivative to 0
      fir_derivative[p] <- 0
    } else {
      # Apply the given formula for the first derivative
      fir_derivative[p] <- (-2*intCol[p-2] - intCol[p-1] + intCol[p+1] + 2*intCol[p+2]) / 10
    }
  }
  
  # Return the first derivative intensities
  return(fir_derivative)
}

calculate_second_derivative <- function(intCol) {
  total_slices <- length(intCol)
  # Initialize second derivative intensities vector
  sec_derivative <- numeric(total_slices)
  
  # Loop through each slice
  for (p in 1:total_slices) {
    if (p <= 2 || p >= total_slices - 1) {
      # For the first two and last two data points, set the second derivative to 0
      sec_derivative[p] <- 0
    } else {
      # Apply the given formula for the second derivative
      sec_derivative[p] <- (2*intCol[p-2] - intCol[p-1] - 2*intCol[p] - intCol[p+1] + 2*intCol[p+2]) / 7
    }
  }
  
  # Return the second derivative intensities
  return(sec_derivative)
}

# Function arPLS_baseline_pos, to do arPLS baseline calculation based on (Baek et al., 2015) Matlab codes 
# Codes are modified so that baseline will always be positive, weight calculation optimized
# lambda: Smoothness penalty parameter and ratio: Convergence tolerance for weight change are suggest values
# of the paper and max_iter = 500 is an emperical iteration maximum limit in case very little optimization could 
# be made even if giving more calcuting time 
arPLS_baseline_pos <- function(y, lambda = 1e5, ratio = 1e-6, max_iter = 500) {
  N <- length(y)
  D <- diff(diag(N), differences = 2)
  H <- lambda * t(D) %*% D
  w <- rep(1, N)
  iter <- 0
  epsilon <- .Machine$double.eps  # Smallest positive number
  
  repeat {
    iter <- iter + 1
    
    W <- diag(w)
    WH <- W + H
    C <- chol(WH)
    
    z <- solve(C, solve(t(C), w * y))
    z[z < 0] <- 0
    
    d <- y - z
    dn <- d[d < 0]
    m <- mean(dn, na.rm = TRUE)
    s <- sd(dn, na.rm = TRUE)
    
    # Handle NA values for mean and standard deviation
    if (is.na(m)) m <- 0
    if (is.na(s) || s < epsilon) s <- epsilon
    
    # Clip the values to avoid overflow in exp
    exp_values <- 2 * (d - (2 * s - m)) / s
    exp_values <- pmin(pmax(exp_values, -500), 500)  # Clipping between -500 and 500
    
    wt <- 1 / (1 + exp(exp_values))
    
    if (any(is.na(wt)) || any(is.infinite(wt))) {
      cat("Encountered NA or Inf in weights. Stopping...\n")
      break
    }
    
    rel_change <- sqrt(sum((w - wt)^2) / sum(w^2))
    cat("Iteration:", iter, "Relative Change:", rel_change, "\n")
    
    if (is.na(rel_change) || rel_change < ratio) {
      cat("Converged after", iter, "iterations.\n")
      break
    }
    
    w <- wt
    
    if (iter >= max_iter) {
      cat("Reached maximum iterations. Stopping...\n")
      break
    }
  }
  
  return(z)
}

# Function check_lengths to do quick check if input lengths are all the same when passing objects to the next function
check_lengths <- function(...) {
  inputs <- list(...)
  lengths <- sapply(inputs, length)
  
  if (length(unique(lengths)) != 1) {
    warning("The lengths of the inputs are not equal")
    return(FALSE)
  } else {
    return(TRUE)
  }
}

# Function convert_List to convert itaList (the T or F list to determine which slice to be kept) of peak region +-180s to the itsList format of +-180s peak 
# 3 minutes before peak start and after peak end are needed to extract information to calculate noise and baseline 
convert_List <- function(List, slices_rt, slices_rt_180s){
  if (!check_lengths(List, slices_rt)) {
    warning("convert_List terminated due to unequal lengths")
    return(NA)
  }
  
  newList <- rep(0, length(slices_rt_180s))
  p_convertion_i <- which(slices_rt_180s == slices_rt[1])
  for (i in 1:length(slices_rt)){
    newList[i + p_convertion_i - 1] <- List[i]
  }
  return(newList)  
}

# Function rv_convert_baseline to get the baseline of peak region +-0s from the baseline of peak region +-180s
rv_convert_baseline <- function(baseline_180s, slices_rt, slices_rt_180s){
  if (!check_lengths(baseline_180s, slices_rt_180s)) {
    warning("rv_convert_baseline terminated due to unequal lengths")
    return(NA)
  }
  
  newbaseline <- rep(0, length(slices_rt))
  p_convertion_i <- which(slices_rt_180s == slices_rt[1])
  for (i in 1:length(slices_rt)){
    newbaseline[i] <- baseline_180s[i + p_convertion_i - 1]
  }
  return(newbaseline)  
}

# Function pick_edges to determine peak edges by back-tracing from peak height to left and right retention sides 
# to make sure both smothered and smoothed slice M0 intensity is >= slice_int_cutoff and also the baseline, be cautious to input the correct baseline region 
pick_edges <- function(intCol, intCol_sm, position, baseline) {
  if (!check_lengths(intCol_sm, baseline)) {
    warning("pick_edges terminated due to unequal lengths")
    return(NA)
  }
  
  total_slices <- length(intCol_sm)
  edgeList <- rep(0, total_slices)
  if (intCol_sm[position] < slice_int_cutoff | intCol[position] < slice_int_cutoff){
    return(edgeList)
  } else {
    edgeList[position] <- 1
  }
  
  # Adjust left edge
  stop_left <- FALSE
  for (p in seq(max(position - 1, 1), 1, by = -1)) {
    if (stop_left) break
    if (intCol_sm[p] <= baseline[p] | intCol_sm[p] < slice_int_cutoff | intCol[p] <= baseline[p] | intCol[p] < slice_int_cutoff) {
      stop_left <- TRUE
    } else {
      edgeList[p] <- 1
    }
  }
  
  # Adjust right edge
  stop_right <- FALSE
  # Adjust right edge
  for (p in seq(min(position + 1, total_slices), total_slices, by = 1)) {
    if (stop_right) break
    if (intCol_sm[p] <= baseline[p] | intCol_sm[p] < slice_int_cutoff | intCol[p] <= baseline[p] | intCol[p] < slice_int_cutoff) {
      stop_right <- TRUE
    } else {
      edgeList[p] <- 1
    }
  }
  
  return(edgeList)
}  

# Function check_continuous_values to check for 5 continuous values > S_N_BL_ratio_cutoff * baseline to 
# know if a slice could be in the left, right, or center part of a peak using also its 1st and 2nd derivative. 
# The info is needed so that only region with no other peaks will be selected to calculate noise.
# 4000 for 1st derivative and 6000 for 2nd derivative is an emperical value and at least work well for Orbitrap
check_continuous_values <- function(start, end, intCol_180s_sm, baseline_180s) {
  fir_der <- calculate_first_derivative(intCol_180s_sm)
  sec_der <- calculate_second_derivative(intCol_180s_sm)
  if (!check_lengths(intCol_180s_sm, baseline_180s)) {
    warning("check_continuous_values terminated due to unequal lengths")
    return(NA)
  }
  
  for (i in seq(start, end - 4)) {  # Subtract 4 to prevent out-of-bounds
    if (all(intCol_180s_sm[i:(i + 4)] > S_N_BL_ratio_cutoff * baseline_180s[i:(i + 4)])) {
      if (fir_der[i+4] > 4000) {
        if (fir_der[i] < fir_der[i+1] & fir_der[i+1] < fir_der[i+2] & fir_der[i+2] < fir_der[i+3] & fir_der[i+3] < fir_der[i+4]) {
          return(TRUE)
        }
      }
      if (all(sec_der[i:(i + 4)] > 4000)) {
        return(TRUE)
      }
      if (fir_der[i] > 4000 & fir_der[i+4] < -4000) {
        if (fir_der[i] > fir_der[i+1] & fir_der[i+1] > fir_der[i+2] & fir_der[i+2] > fir_der[i+3] & fir_der[i+3] > fir_der[i+4]) {
          if (any(sec_der[i:(i + 4)] < -6000)) {
            return(TRUE)
          } 
        }
      }
      if (all(sec_der[i:(i + 4)] < -4000)) {
        return(TRUE)
      }
      if (fir_der[i] < -4000) {
        if (fir_der[i] < fir_der[i+1] & fir_der[i+1] < fir_der[i+2] & fir_der[i+2] < fir_der[i+3] & fir_der[i+3] < fir_der[i+4]) {
          return(TRUE)
        } 
      }
    }
  }
  
  return(FALSE)
}


# Function estimate_noise to use RMS based method to calculate noise by checking the 180 seconds region before the metabolite peak, 
# and find a region of 30 seconds without other peaks to calculate noise. If there is always another metabolite peak, a 10 second window is used instead
# to calculate noise. In rare ocasions when a 10 second window is also unavailable, noise will be NA
estimate_noise <- function(edgeList_180s, intCol_180s, intCol_180s_sm, slices_rt_180s, baseline_180s) {
  if (!check_lengths(edgeList_180s, intCol_180s, intCol_180s_sm, slices_rt_180s, baseline_180s)) {
    warning("estimate_noise terminated due to unequal lengths")
    return(NA)
  }
  
  if (all(edgeList_180s == 0)) {
    return(NA)
  }
  
  p_right <- which(edgeList_180s == 1)[1] - 1
  if (is.na(p_right)) {
    return(NA)
  }
  if (slices_rt_180s[p_right] - 30 < 0) {
    warning(paste(compound, ": region selected to calculate noise is too close to 0 min thus become less than 30s, 10s is used instead"))
    return(estimate_noise_10s(edgeList_180s, intCol_180s, intCol_180s_sm, slices_rt_180s, baseline_180s))
  }
  
  rtime_left <- max(0, slices_rt_180s[p_right] - 30)
  p_left <- which(slices_rt_180s >= rtime_left)[1]
  
  repeat {
    p_avoid <- check_continuous_values(p_left, p_right, intCol_180s_sm, baseline_180s)
    
    if (!(p_avoid)) {
      break
    } else {
      p_right <- p_avoid - 1
      if (p_right <= 0) {
        warning(paste(compound, ": shifted region selected to calculate noise is too close to 0 min thus become less than 30s, 10s is used instead"))
        return(estimate_noise_10s(edgeList_180s, intCol_180s, intCol_180s_sm, slices_rt_180s, baseline_180s))
      }
      
      if (slices_rt_180s[p_right] - 30 < 0) {
        warning(paste(compound, ": shifted region selected to calculate noise is too close to 0 min thus become less than 30s, 10s is used instead"))
        return(estimate_noise_10s(edgeList_180s, intCol_180s, intCol_180s_sm, slices_rt_180s, baseline_180s))
      }
      
      rtime_left <- max(0, slices_rt_180s[p_right] - 30)
      p_left <- which(slices_rt_180s >= rtime_left)[1]
      
      if (is.na(p_left) || p_left <= 0) {
        warning(paste(compound, ": shifted region selected to calculate noise is too close to 0 min thus become less than 30s, 10s is used instead"))
        return(estimate_noise_10s(edgeList_180s, intCol_180s, intCol_180s_sm, slices_rt_180s, baseline_180s))
      }
    }
  }
  
  difference <- intCol_180s[p_left:p_right] - baseline_180s[p_left:p_right]
  cat("\nestimating noises using 30s region.")
  cat("\np_left is", p_left, ". ")
  cat("\np_right is", p_right, ". \n")
  n <- length(difference)
  # n-2 is used by correcting the RMS method bias when introducing baseline
  rms_noise <- sqrt(sum(difference^2) / (n - 2))
  return(rms_noise)
}

estimate_noise_10s <- function(edgeList_180s, intCol_180s, intCol_180s_sm, slices_rt_180s, baseline_180s) {
  if (!check_lengths(edgeList_180s, intCol_180s, intCol_180s_sm, slices_rt_180s, baseline_180s)) {
    warning("estimate_noise terminated due to unequal lengths")
    return(NA)
  }
  
  if (all(edgeList_180s == 0)) {
    return(NA)
  }
  
  p_right <- which(edgeList_180s == 1)[1] - 1
  if (is.na(p_right)) {
    return(NA)
  }
  if (slices_rt_180s[p_right] - 10 < 0) {
    warning(paste(compound, ": region selected to calculate noise is too close to 0 min thus become less than 10s, terminated"))
    return(NA)
  }
  
  rtime_left <- max(0, slices_rt_180s[p_right] - 10)
  p_left <- which(slices_rt_180s >= rtime_left)[1]
  
  repeat {
    p_avoid <- check_continuous_values(p_left, p_right, intCol_180s_sm, baseline_180s)
    
    if (!(p_avoid)) {
      break
    } else {
      p_right <- p_avoid - 1
      if (p_right <= 0) {
        warning(paste(compound, ": no 10s region contains no peak to be used for noise calcualtion, terminated"))
        return(NA)
      }
      
      if (slices_rt_180s[p_right] - 10 < 0) {
        warning(paste(compound, ": no 10s region contains no peak to be used for noise calcualtion, terminated"))
        return(NA)
      }
      
      rtime_left <- max(0, slices_rt_180s[p_right] - 10)
      p_left <- which(slices_rt_180s >= rtime_left)[1]
      
      if (is.na(p_left) || p_left <= 0) {
        warning(paste(compound, ": no 10s region contains no peak to be used for noise calcualtion, terminated"))
        return(NA)
      }
    }
  }
  
  difference <- intCol_180s[p_left:p_right] - baseline_180s[p_left:p_right]
  cat("\nestimating noises using 10s region.")
  cat("\np_left is", p_left, ". ")
  cat("\np_right is", p_right, ". \n")
  n <- length(difference)
  # n-2 is used by correcting the RMS method bias when introdcuing baseline
  rms_noise <- sqrt(sum(difference^2) / (n - 2))
  return(rms_noise)
}

# Function check_elements_800, check if all elements > correlation_cutoff
check_elements_800 <- function(correlation_half_part, f800, l800) {
  return(all(correlation_half_part[f800:l800] > correlation_cutoff))
}

# Function check_peak to check if a peak is indeed a peak by checking for 5 continuous values > S_N_BL_ratio_cutoff * baseline, 
# and 1st derivative decrease from 4000 to -4000, and 2nd derivative any < -6000  
check_peak <- function(start, end, intCol_sm, baseline) {
  fir_der <- calculate_first_derivative(intCol_sm)
  sec_der <- calculate_second_derivative(intCol_sm)
  if (!check_lengths(intCol_sm, baseline)) {
    warning("check_peak terminated due to unequal lengths")
    return(NA)
  }
  
  for (i in seq(start, end - 4)) {  # Subtract 4 to prevent out-of-bounds
    if (all(intCol_sm[i:(i + 4)] > S_N_BL_ratio_cutoff * baseline[i:(i + 4)])) {
      if (fir_der[i] > 4000 & fir_der[i+4] < -4000) {
        if (fir_der[i] > fir_der[i+1] & fir_der[i+1] > fir_der[i+2] & fir_der[i+2] > fir_der[i+3] & fir_der[i+3] > fir_der[i+4]) {
          if (any(sec_der[i:(i + 4)] < -6000)) {
            return(TRUE)
          } 
        }
      }
    }
  }
  return(FALSE)
}

# Function localMaxima to find local maxima using the provided code, 
# from https://stackoverflow.com/questions/6836409/finding-local-maxima-and-minima
localMaxima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  y <- diff(c(-.Machine$integer.max, x)) > 0L
  rle_y <- rle(y)$lengths
  y <- cumsum(rle_y)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  y
}

# Function find_highest_local_max_in_regions to find the highest local maximum within specified regions that are peaks
find_highest_local_max_in_regions <- function(intCol_sm, baseline, position, left_edge_p, right_edge_p) {
  if (!check_lengths(intCol_sm, baseline)) {
    warning("find_highest_local_max_in_regions")
    return(NA)
  }
  total_slices <- length(intCol_sm)
  # Find all local maxima
  local_maxima_positions <- localMaxima(intCol_sm)
  
  # Define the regions
  region1 <- which(local_maxima_positions > left_edge_p & local_maxima_positions < position - 3)
  region2 <- which(local_maxima_positions > position + 3 & local_maxima_positions < right_edge_p)
  
  # Combine the regions
  valid_positions <- c(local_maxima_positions[region1], local_maxima_positions[region2])
  
  # Check if there are valid positions
  if (length(valid_positions) == 0) {
    return(NA)
  } else {
    final_valid_positions <- c()
    # Find the position of the highest value among the valid positions
    for (p in valid_positions) {
      start <- p - 2
      end <- p + 2
      if (start < 1){
        start <- 1
        end <- 5
      }
      if (end > total_slices){
        start <- total_slices - 4
        end <- total_slices
      }
      if (check_peak(start, end, intCol_sm, baseline)){
        final_valid_positions <- c(final_valid_positions, p)
      }
    }
  }
  if (length(final_valid_positions) == 0) {
    return(NA)
  } else {
    max_value_position <- final_valid_positions[which.max(intCol_sm[final_valid_positions])]
    return(max_value_position)
  }
}

# Function pick_slices that pick high quality slices for M0M1M2 ratio calculations from the initial determine peak slices
# based on slice_int_cutoff, S_N_BL_ratio_cutoff, and correlation_cutoff
pick_slices <- function(edgeList, intCol, intCol_sm, position, baseline, noise, integrated_to_area_Mprevious, intCOlM0_sm, baselineM0) {
  if (!check_lengths(edgeList, intCol, intCol_sm, baseline, integrated_to_area_Mprevious, intCOlM0_sm)) {
    warning("pick_slices terminated due to unequal lengths")
    return(NA)
  }
  if (is.na(noise)) {
    warning("noise is NA and treated as 0")
    noise <- 0
  }
  
  total_slices <- length(intCol)
  itaList <- rep(0, total_slices)
  
  if (all(edgeList == 0)) {
    itaList <- rep(0, total_slices)
    return(itaList)
  }
  
  left_edge_p <- which(edgeList == 1)[1]
  right_edge_p <- tail(which(edgeList == 1), n = 1)
  
  stop_left <- FALSE
  stop_right <- FALSE
  
  if (intCol[position] < slice_int_cutoff | intCol_sm[position] < slice_int_cutoff | intCol[position] < correlation_cutoff * noise | intCol[position] < correlation_cutoff * baseline[position] | intCol_sm[position] < correlation_cutoff * noise | intCol_sm[position] < correlation_cutoff * baseline[position]) {
    stop_left <- TRUE
    stop_right <- TRUE
    itaList <- rep(0, total_slices)
    warningmessage <<- paste(warningmessage, "Peak intensity is too low, or baseline and noise is too high")
    return(itaList)
  } else {
    itaList[position] <- 1
  }
  
  # Adjust left edge
  for (p in seq(max(position - 1, left_edge_p), left_edge_p, by = -1)) {
    if (stop_left) break
    if (intCol[p] < slice_int_cutoff | intCol_sm[p] < slice_int_cutoff | intCol[p] < S_N_BL_ratio_cutoff * noise | intCol[p] < S_N_BL_ratio_cutoff * baseline[p] | intCol_sm[p] < S_N_BL_ratio_cutoff * noise | intCol_sm[p] < S_N_BL_ratio_cutoff * baseline[p]) {
      stop_left <- TRUE
    } else {
      itaList[p] <- 1
    }
  }
  
  # Adjust right edge
  for (p in seq(min(position + 1, right_edge_p), right_edge_p, by = 1)) {
    if (stop_right) break
    if (intCol[p] < slice_int_cutoff | intCol_sm[p] < slice_int_cutoff | intCol[p] < S_N_BL_ratio_cutoff * noise | intCol[p] < S_N_BL_ratio_cutoff * baseline[p] | intCol_sm[p] < S_N_BL_ratio_cutoff * noise | intCol_sm[p] < S_N_BL_ratio_cutoff * baseline[p]) {
      stop_right <- TRUE
    } else {
      itaList[p] <- 1
    }
  }
  
  # Mn shall width shall never exceeds M0
  if (length(integrated_to_area_Mprevious) > 0) {
    zero_indices <- which(integrated_to_area_Mprevious == 0)
    itaList[zero_indices] <- 0
  }
  
  if(!any(unlist(itaList) == 1)){
    itaList <- rep(0, total_slices)
    return(itaList)
  }
  
  left_edge_p <- which(edgeList == 1)[1]
  right_edge_p <- tail(which(edgeList == 1), n = 1)
  
  data <- intCol_sm[left_edge_p:right_edge_p] - baseline[left_edge_p:right_edge_p]
  reference_shape <- intCOlM0_sm[left_edge_p:right_edge_p] - baselineM0[left_edge_p:right_edge_p]
  
  score_whole_peak <- max(cor(data, reference_shape, method = "pearson"), cor(data, reference_shape, method = "spearman"))
  
  cat("\npicking slices of compound", compound, ". ")
  cat("\ndata is", data, ". ")
  cat("\nreference_shape is", reference_shape, ". ")
  cat("\nwhole peak score is:", score_whole_peak, ". ")
  
  left_slice_p <- which(itaList == 1)[1]
  right_slice_p <- tail(which(itaList == 1), n = 1)
  
  data <- intCol_sm[left_slice_p:right_slice_p] - baseline[left_slice_p:right_slice_p]
  reference_shape <- intCOlM0_sm[left_slice_p:right_slice_p] - baselineM0[left_slice_p:right_slice_p]
  
  score <- max(cor(data, reference_shape, method = "pearson"), cor(data, reference_shape, method = "spearman"))
  
  cat("\ndata is", data, ". ")
  cat("\nreference_shape is", reference_shape, ". ")
  cat("\nscore is:", score, ". ")
  
  if(is.na(score)){
    itaList <- rep(0, total_slices)
    return(itaList)
  }
  
  if(score > correlation_cutoff){
    return(itaList)
  } 
  
  itaList <- rep(0, total_slices)
  return(itaList)
}




################################


