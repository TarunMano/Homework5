
setClass(
  Class = "sparse_numeric",
  slots = c(
    value = "numeric", 
    pos = "integer", 
    length = "integer"))

setValidity("sparse_numeric", function(object) {
  errs = character(0)
  
  if (!is.numeric(object@value))  errs = c(errs, "'value' must be numeric")
  if (!is.integer(object@pos))    errs = c(errs, "'pos' must be integer")
  if (!is.integer(object@length) || length(object@length) != 1L)
    errs = c(errs, "'length' must be a single integer")
  
  n = if (length(object@length)) object@length[1L] else NA_integer_
  
  if (length(object@value) != length(object@pos))
    errs = c(errs, "'value' and 'pos' lengths must match")
  
  if (length(object@pos)) {
    if (any(is.na(object@pos)))              errs = c(errs, "'pos' has NA")
    if (any(object@pos < 1L | object@pos > n)) errs = c(errs, "'pos' outside 1..length")
    if (anyDuplicated(object@pos))           errs = c(errs, "'pos' has duplicates")
    if (!all(diff(object@pos) > 0L))         errs = c(errs, "'pos' must be strictly ascending")
  }
  
  if (!is.na(n) && n < 0L) errs = c(errs, "'length' must be >= 0")
  
  if (length(errs)) errs else TRUE
})

setAs("numeric", "sparse_numeric", function(from) {
  if (!length(from)) {
    return(new("sparse_numeric", value = numeric(), pos = integer(), length = 0L))
  }
  nz  = which(from != 0)
  new("sparse_numeric",
      value  = as.numeric(from[nz]),
      pos    = as.integer(nz),
      length = as.integer(length(from)))
})

setAs("sparse_numeric", "numeric", function(from) {
  out = numeric(from@length)
  if (length(from@pos)) out[from@pos] = from@value
  out
})

setMethod("show", "sparse_numeric", function(object) {
  nnz = length(object@pos)
  cat("sparse_numeric (length =", object@length, ", nnz =", nnz, ")\n")
  if (nnz) {
    k  = min(nnz, 8L)
    df = data.frame(pos = object@pos[seq_len(k)],
                     value = object@value[seq_len(k)])
    print(df, row.names = FALSE)
    if (nnz > k) cat("  ...", nnz - k, "more non-zeros\n")
  } else {
    cat("  <all zeros>\n")
  }
})

setGeneric("sparse_add",function(x, y, ...) standardGeneric("sparse_add"))
setGeneric("sparse_sub",function(x, y, ...) standardGeneric("sparse_sub"))
setGeneric("sparse_mult",function(x, y, ...) standardGeneric("sparse_mult"))
setGeneric("sparse_crossprod",function(x, y, ...) standardGeneric("sparse_crossprod"))

.same_len = function(x, y) {
  if (x@length != y@length) stop("Sparse vectors must have the same 'length'.")
}


setMethod("sparse_add", c("sparse_numeric", "sparse_numeric"), function(x, y) {
  .same_len(x, y)
  if (!length(x@pos)) return(x = new("sparse_numeric", value = y@value, pos = y@pos, length = y@length))
  if (!length(y@pos)) return(x)
  
  all_pos = sort(unique(c(x@pos, y@pos)))
  
  vx = numeric(length(all_pos))
  vy = numeric(length(all_pos))
  
  vx[match(x@pos, all_pos)] = x@value
  vy[match(y@pos, all_pos)] = y@value
  
  v  = vx + vy
  keep = v != 0
  new("sparse_numeric",
      value = v[keep],
      pos   = as.integer(all_pos[keep]),
      length = x@length)
})

setMethod("sparse_sub", c("sparse_numeric", "sparse_numeric"), function(x, y) {
  .same_len(x, y)
  if (!length(y@pos)) return(x)
  if (!length(x@pos)) return(new("sparse_numeric",
                                 value = -y@value, pos = y@pos, length = y@length))
  
  all_pos = sort(unique(c(x@pos, y@pos)))
  
  vx = numeric(length(all_pos))
  vy = numeric(length(all_pos))
  
  vx[match(x@pos, all_pos)] = x@value
  vy[match(y@pos, all_pos)] = y@value
  
  v  = vx - vy
  keep = v != 0
  new("sparse_numeric",
      value = v[keep],
      pos   = as.integer(all_pos[keep]),
      length = x@length)
})

setMethod("sparse_mult", c("sparse_numeric", "sparse_numeric"), function(x, y) {
  .same_len(x, y)
  if (!length(x@pos) || !length(y@pos))
    return(new("sparse_numeric", value = numeric(), pos = integer(), length = x@length))
  
  both = intersect(x@pos, y@pos)
  if (!length(both))
    return(new("sparse_numeric", value = numeric(), pos = integer(), length = x@length))
  
  ix = match(both, x@pos)
  iy = match(both, y@pos)
  v  = x@value[ix] * y@value[iy]
  
  keep = v != 0
  new("sparse_numeric",
      value = v[keep],
      pos   = as.integer(both[keep]),
      length = x@length)
})

setMethod("sparse_crossprod", c("sparse_numeric", "sparse_numeric"), function(x, y) {
  .same_len(x, y)
  if (!length(x@pos) || !length(y@pos)) return(0.0)
  both = intersect(x@pos, y@pos)
  if (!length(both)) return(0.0)
  ix = match(both, x@pos)
  iy = match(both, y@pos)
  sum(x@value[ix] * y@value[iy])
})

setMethod("+", c("sparse_numeric", "sparse_numeric"),
          function(e1, e2) sparse_add(e1, e2))
setMethod("-", c("sparse_numeric", "sparse_numeric"),
          function(e1, e2) sparse_sub(e1, e2))
setMethod("*", c("sparse_numeric", "sparse_numeric"),
          function(e1, e2) sparse_mult(e1, e2))

setMethod("plot", c("sparse_numeric", "sparse_numeric"), function(x, y, ...) {
  .same_len(x, y)
  
  xmin = 1; xmax = max(1L, x@length)
  vals = c(0, x@value, y@value)
  ymin = min(vals); ymax = max(vals)
  if (ymin == ymax) { ymin = ymin - 1; ymax = ymax + 1 }
  
  plot(NA, xlim = c(xmin, xmax), ylim = c(ymin, ymax),
       xlab = "Position", ylab = "Value",
       main = "Sparse non-zeros", ...)
  
  if (length(x@pos)) points(x@pos, x@value, pch = 16)
  if (length(y@pos)) points(y@pos, y@value, pch = 1)
  
  ov = intersect(x@pos, y@pos)
  if (length(ov)) {
    xv = x@value[match(ov, x@pos)]
    yv = y@value[match(ov, y@pos)]
    segments(ov, xv, ov, yv, lty = 2)
  }
  
  legend("topleft", bty = "n",
         legend = c("x (filled)", "y (open)", "overlap"),
         pch = c(16, 1, NA), lty = c(NA, NA, 2))
})

setMethod("length", "sparse_numeric", function(x) x@length)
