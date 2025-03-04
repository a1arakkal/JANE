
trace_plot <- function(model_object){
  
  plot_data <- model_object$termination_metric
  
  termination_rule_labels <- data.frame(termination_rule = c("ARI",
                                                             "NMI",
                                                             "CER",
                                                             "prob_mat",
                                                             "Q"),
                                        label = c("ARI current vs. previous iteration",
                                                  "NMI current vs. previous iteration",
                                                  "CER current vs. previous iteration",
                                                  "Cluster assignment probability matrix",
                                                  "Abs diff. Q current vs. previous iteration"))
  
  opar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar))
  
  if(model_object$termination_rule == "prob_mat"){
    
    graphics::layout(matrix(c(1,1,1, 
                              2,2,2, 
                              3,3, 
                              4,4,
                              5,5), 
                            nrow = 2, byrow = TRUE), heights = c(1, 0.7))
    
    temp_line_lab <- strsplit(colnames(plot_data)[4], "_*", fixed = T)[[1]]
    has_prob_W <- length(grep(pattern = "prob_mat_W", x = colnames(plot_data)))>0
    
    graphics::par(mar = c(5.1, 5.1, 1.1, 2.1))
    plot(plot_data[,4], type = "l",
         xlab = "Iteration",
         ylab = ifelse(!has_prob_W,
                       paste0("Abs diff. cluster assignment prob. matrix\n current vs. previous iteration",
                              ifelse(temp_line_lab[2] == "1", " (max)",
                                     paste0(" (quant. = ", as.numeric(temp_line_lab[2]), ")"))),
                       paste0("Abs diff. current vs. previous iteration\n ",
                              ifelse(temp_line_lab[2] == "1", " (max)",
                                     paste0(" (quant. = ", as.numeric(temp_line_lab[2]), ")")))))
    if(has_prob_W){
      graphics::lines(plot_data[,5], col = "blue", type = "l")
      graphics::legend("topright", 
                       legend=c("Actor cluster assignment prob. matrix", "Edge weight cluster assignment prob. matrix"),
                       col=c("black", "blue"), lty = 1, cex=0.8)
    }
    
    graphics::par(mar = c(5.1, 4.1, 1.1, 2.1))
    finite_vals <- which(is.finite(plot_data[, grep("abs_diff_MA_U", colnames(plot_data))]))
    
    if(length(finite_vals)>1){
      plot(y = plot_data[finite_vals, grep("abs_diff_MA_U", colnames(plot_data))], 
           x = finite_vals,
           type = "l",
           xlab = "Iteration",
           ylab = "Abs diff. cumulative average",
           col = "red", 
           ylim = if(!has_prob_W){
             c(min(plot_data[finite_vals, grep("abs_diff_MA_U", colnames(plot_data))-1],
                   plot_data[finite_vals,grep("abs_diff_MA_U", colnames(plot_data))]), 
               max(plot_data[finite_vals, grep("abs_diff_MA_U", colnames(plot_data))-1], 
                   plot_data[finite_vals,grep("abs_diff_MA_U", colnames(plot_data))]))
           } else{
             c(min(plot_data[finite_vals, grep("abs_diff_MA_U", colnames(plot_data))-1],
                   plot_data[finite_vals, grep("abs_diff_MA_U", colnames(plot_data))-2],
                   plot_data[finite_vals,grep("abs_diff_MA_U", colnames(plot_data))]), 
               max(plot_data[finite_vals, grep("abs_diff_MA_U", colnames(plot_data))-1], 
                   plot_data[finite_vals, grep("abs_diff_MA_U", colnames(plot_data))-2],
                   plot_data[finite_vals,grep("abs_diff_MA_U", colnames(plot_data))]))
           }
      )
      if(!has_prob_W){
        graphics::lines(plot_data[, grep("abs_diff_MA_U", colnames(plot_data))-1], col = "blue", type = "l")
        graphics::legend("topright", 
                         legend=c("U", "Cluster assignment prob. matrix"),
                         col=c("red", "blue"), lty = 1, cex=0.8)
      } else {
        graphics::lines(plot_data[, grep("abs_diff_MA_U", colnames(plot_data))-2], col = "black", type = "l")
        graphics::lines(plot_data[, grep("abs_diff_MA_U", colnames(plot_data))-1], col = "blue", type = "l")
        graphics::legend("topright", 
                         legend=c("U", "Actor cluster assignment prob. matrix", "Edge weight cluster assignment prob. matrix"),
                         col=c("red", "black" ,"blue"), lty = 1, cex=0.8)
      }
    } else {
      plot(0,
           xlim = c(0, nrow(plot_data)),
           yaxt = 'n', pch = '',
           xlab = "Iteration",
           ylab = "Abs diff. cumulative average")
      if(!has_prob_W){
        graphics::legend("topright", 
                         legend=c("U", "Cluster assignment prob. matrix"),
                         col=c("red", "blue"), lty = 1, cex=0.8)
      }else {
        graphics::legend("topright", 
                         legend=c("U", "Actor cluster assignment prob. matrix", "Edge weight cluster assignment prob. matrix"),
                         col=c("red", "black" ,"blue"), lty = 1, cex=0.8)
      }
    }
    
    graphics::par(mar = c(5.1, 4.1, 1.1, 2.1))
    plot(plot_data[,"ARI"], type = "l",
         xlab = "Iteration",
         ylab = "ARI current vs. previous iteration",
         col = "blue")
    
    graphics::par(mar = c(5.1, 4.1, 1.1, 2.1))
    plot(plot_data[,"NMI"], type = "l",
         xlab = "Iteration",
         ylab = "NMI current vs. previous iteration",
         col = "darkgreen")
    
    graphics::par(mar = c(5.1, 4.1, 1.1, 2.1))
    plot(plot_data[,"CER"], type = "l",
         xlab = "Iteration",
         ylab = "CER current vs. previous iteration",
         col = "red")
    
  } else if (model_object$termination_rule == "Q") {
    
    graphics::layout(matrix(c(1,1,1, 
                              2,2,2, 
                              3,3, 
                              4,4,
                              5,5), 
                            nrow = 2, byrow = TRUE), heights = c(1, 0.7))
    
    graphics::par(mar = c(5.1, 4.1, 1.1, 2.1))
    plot(plot_data[,4], type = "l",
         xlab = "Iteration",
         ylab = "Abs diff. Q current vs. previous iteration")
    
    graphics::par(mar = c(5.1, 4.1, 1.1, 2.1))
    finite_vals <- which(is.finite(plot_data[,5]))
    
    if (length(finite_vals)>1){
      plot(y = plot_data[finite_vals,5], 
           x = finite_vals,
           type = "l",
           xlab = "Iteration",
           ylab = "Abs diff. Q cumulative average")
    } else {
      plot(0,
           xlim = c(0, nrow(plot_data)),
           yaxt = 'n', pch = '',
           xlab = "Iteration",
           ylab = "Abs diff. Q cumulative average")
    }
    
    graphics::par(mar = c(5.1, 4.1, 1.1, 2.1))
    plot(plot_data[,"ARI"], type = "l",
         xlab = "Iteration",
         ylab = "ARI current vs. previous iteration",
         col = "blue")
    
    graphics::par(mar = c(5.1, 4.1, 1.1, 2.1))
    plot(plot_data[,"NMI"], type = "l",
         xlab = "Iteration",
         ylab = "NMI current vs. previous iteration",
         col = "darkgreen")
    
    graphics::par(mar = c(5.1, 4.1, 1.1, 2.1))
    plot(plot_data[,"CER"], type = "l",
         xlab = "Iteration",
         ylab = "CER current vs. previous iteration",
         col = "red")
    
  } else {
    
    graphics::layout(matrix(c(1,1, 
                              2,2,
                              3,3, 
                              4,4), 
                            nrow = 2, byrow = TRUE), heights = c(1, 1))
    
    graphics::par(mar = c(5.1, 4.1, 1.1, 2.1))
    finite_vals <- which(is.finite(plot_data[,6]))
    
    if(length(finite_vals)>1){
      plot(y = plot_data[finite_vals,6], 
           x = finite_vals,
           type = "l",
           xlab = "Iteration",
           ylab = "Abs diff. cumulative average",
           col = "red", 
           ylim = c(min(plot_data[finite_vals,5], plot_data[finite_vals,6]), 
                    max(plot_data[finite_vals,5], plot_data[finite_vals,6])))
      graphics::lines(plot_data[,5], col = "blue", type = "l")
      graphics::legend("topright", legend=c("U", model_object$termination_rule),
                       col=c("red", "blue"), lty = 1, cex=0.8)
    } else {
      plot(0,
           xlim = c(0, nrow(plot_data)),
           yaxt = 'n', pch = '',
           xlab = "Iteration",
           ylab = "Abs diff. cumulative average")
      graphics::legend("topright", legend=c("U", model_object$termination_rule),
                       col=c("red", "blue"), lty = 1, cex=0.8)
    }
    
    graphics::par(mar = c(5.1, 4.1, 1.1, 2.1))
    plot(plot_data[,"ARI"], type = "l",
         xlab = "Iteration",
         ylab = "ARI current vs. previous iteration",
         col = "blue")
    
    graphics::par(mar = c(5.1, 4.1, 1.1, 2.1))
    plot(plot_data[, "NMI"], type = "l",
         xlab = "Iteration",
         ylab = "NMI current vs. previous iteration",
         col = "darkgreen")
    
    graphics::par(mar = c(5.1, 4.1, 1.1, 2.1))
    plot(plot_data[, "CER"], type = "l",
         xlab = "Iteration",
         ylab = "CER current vs. previous iteration",
         col = "red")
    
  }
  
}
