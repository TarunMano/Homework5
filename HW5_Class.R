## HW5 Class/Methods

setClass(
  Class = "sparse_numeric",
  slots = c(
    value  = "numeric",
    pos    = "integer",
    length = "integer"
  )
)

## -------------------------------
## Validity method
## -------------------------------
setValidity("sparse_numeric", function(object) {
  errs <- character()
  
  ## slot types (S4 already enforces, but we double-check with messages)
  if (!is.numeric(object@value)) errs <- c(errs, "'value' must be numeric.")
  if (!is.integer(object@pos))   errs <- c(errs,   "'pos' must be integer.")
  if (!is.integer(object@length) || length(object@length) != 1L)
    errs <- c(errs, "'length' must be a single integer.")
  
  ## structural checks
  if (length(object@value) != length(object@pos)) {
    errs <- c(errs, "'value' and 'pos' must have the same length.")
  }
  
  n <- if (length(object@length)) object@length[1L] else NA_integer_
  
  if (length(object@pos)) {
    if (any(is.na(object@pos))) errs <- c(errs, "'pos' contains NA.")
    if (any(object@pos < 1L | object@pos > n))
      errs <- c(errs, "'pos' must lie in 1..length.")
    if (anyDuplicated(object@pos))
      errs <- c(errs, "'pos' must not contain duplicates.")
    ## optional but helpful: enforce sorted positions
    if (!isTRUE(all(diff(object@pos) > 0L)))
      errs <- c(errs, "'pos' must be strictly ascending.")
  }
  
  ## non-negative length
  if (!is.na(n) && n < 0L) errs <- c(errs, "'length' must be >= 0.")
  
  if (length(errs)) errs else TRUE
})

## -------------------------------
## Helper constructors & utilities
## -------------------------------

## Build a sparse_numeric safely (sorts, drops ~zero, validates)
.sparse_build <- function(value, pos, length, tol = sqrt(.Machine$double.eps)) {
  stopifnot(is.numeric(value), is.integer(pos), is.integer(length), length(length) == 1L)
  
  if (length(value) != length(pos)) stop("value and pos must have equal length.")
  if (length(value)) {
    keep <- abs(value) > tol
    value <- value[keep]
    pos   <- pos[keep]
    if (length(pos)) {
      o <- order(pos)
      pos <- pos[o]
      value <- value[o]
    }
  }
  new("sparse_numeric", value = value, pos = pos, length = length)
}

## Check same overall length
.check_same_len <- function(x, y) {
  if (x@length != y@length) stop("Sparse vectors must have the same 'length'.")
}

## Merge-by-index utilities (no dense expansion)
## Returns union/intersection of indices with left/right matches
.align_union <- function(x, y) {
  idx <- sort(unique(c(x@pos, y@pos)))
  ix  <- match(idx, x@pos, nomatch = 0L)
  iy  <- match(idx, y@pos, nomatch = 0L)
  list(idx = idx, ix = ix, iy = iy)
}

.align_intersect <- function(x, y) {
  idx <- intersect(x@pos, y@pos)
  if (length(idx)) idx <- sort(idx)
  ix <- match(idx, x@pos)
  iy <- match(idx, y@pos)
  list(idx = idx, ix = ix, iy = iy)
}

## -------------------------------
## Coercions
## -------------------------------

## numeric -> sparse_numeric
setAs("numeric", "sparse_numeric", function(from) {
  if (!length(from)) {
    return(.sparse_build(numeric(), integer(), 0L))
  }
  pos <- which(from != 0)
  .sparse_build(
    value  = as.numeric(from[pos]),
    pos    = as.integer(pos),
    length = as.integer(length(from))
  )
})

## sparse_numeric -> numeric
setAs("sparse_numeric", "numeric", function(from) {
  out <- numeric(from@length)
  if (length(from@pos)) out[from@pos] <- from@value
  out
})

## -------------------------------
## Show method
## -------------------------------
setMethod("show", "sparse_numeric", function(object) {
  nnz <- length(object@pos)
  cat("sparse_numeric (length =", object@length, ", nnz =", nnz, ")\n", sep = " ")
  if (nnz == 0L) {
    cat("  <all zeros>\n")
  } else {
    max_print <- 8L
    k <- min(nnz, max_print)
    df <- data.frame(pos = object@pos[seq_len(k)], value = object@value[seq_len(k)])
    print(df, row.names = FALSE)
    if (nnz > max_print) cat("  ...", nnz - max_print, "more non-zeros\n")
  }
})

## -------------------------------
## Generics
## -------------------------------
setGeneric("sparse_add",   function(x, y, ...) standardGeneric("sparse_add"))
setGeneric("sparse_sub",   function(x, y, ...) standardGeneric("sparse_sub"))
setGeneric("sparse_mult",  function(x, y, ...) standardGeneric("sparse_mult"))
setGeneric("sparse_crossprod",
           function(x, y, ...) standardGeneric("sparse_crossprod"))

## -------------------------------
## Methods: arithmetic (sparse-only)
## -------------------------------

## Addition: union of positions, sum values, drop ~zeros
setMethod("sparse_add",
          signature(x = "sparse_numeric", y = "sparse_numeric"),
          function(x, y, tol = sqrt(.Machine$double.eps)) {
            .check_same_len(x, y)
            if (length(x@pos) == 0L) return(.sparse_build(y@value, y@pos, y@length))
            if (length(y@pos) == 0L) return(.sparse_build(x@value, x@pos, x@length))
            
            a <- .align_union(x, y)
            
            vx <- numeric(length(a$idx))
            vy <- numeric(length(a$idx))
            selx <- a$ix > 0L
            sely <- a$iy > 0L
            if (any(selx)) vx[selx] <- x@value[a$ix[selx]]
            if (any(sely)) vy[sely] <- y@value[a$iy[sely]]
            
            v <- vx + vy
            keep <- abs(v) > tol
            .sparse_build(v[keep], as.integer(a$idx[keep]), x@length)
          })

## Subtraction: union of positions, subtract values, drop ~zeros
setMethod("sparse_sub",
          signature(x = "sparse_numeric", y = "sparse_numeric"),
          function(x, y, tol = sqrt(.Machine$double.eps)) {
            .check_same_len(x, y)
            if (length(y@pos) == 0L) return(.sparse_build(x@value, x@pos, x@length))
            if (length(x@pos) == 0L) return(.sparse_build(-y@value, y@pos, y@length))
            
            a <- .align_union(x, y)
            
            vx <- numeric(length(a$idx))
            vy <- numeric(length(a$idx))
            selx <- a$ix > 0L
            sely <- a$iy > 0L
            if (any(selx)) vx[selx] <- x@value[a$ix[selx]]
            if (any(sely)) vy[sely] <- y@value[a$iy[sely]]
            
            v <- vx - vy
            keep <- abs(v) > tol
            .sparse_build(v[keep], as.integer(a$idx[keep]), x@length)
          })

## Elementwise multiplication: intersection of positions
setMethod("sparse_mult",
          signature(x = "sparse_numeric", y = "sparse_numeric"),
          function(x, y, tol = sqrt(.Machine$double.eps)) {
            .check_same_len(x, y)
            if (length(x@pos) == 0L || length(y@pos) == 0L) {
              return(.sparse_build(numeric(), integer(), x@length))
            }
            a <- .align_intersect(x, y)
            if (!length(a$idx)) return(.sparse_build(numeric(), integer(), x@length))
            v <- x@value[a$ix] * y@value[a$iy]
            keep <- abs(v) > tol
            .sparse_build(v[keep], as.integer(a$idx[keep]), x@length)
          })

## Cross product (dot product): numeric scalar
setMethod("sparse_crossprod",
          signature(x = "sparse_numeric", y = "sparse_numeric"),
          function(x, y) {
            .check_same_len(x, y)
            if (length(x@pos) == 0L || length(y@pos) == 0L) return(0.0)
            a <- .align_intersect(x, y)
            if (!length(a$idx)) return(0.0)
            sum(x@value[a$ix] * y@value[a$iy])
          })

## -------------------------------
## Operators +, -, *
## -------------------------------
setMethod("+", signature(e1 = "sparse_numeric", e2 = "sparse_numeric"),
          function(e1, e2) sparse_add(e1, e2))

setMethod("-", signature(e1 = "sparse_numeric", e2 = "sparse_numeric"),
          function(e1, e2) sparse_sub(e1, e2))

setMethod("*", signature(e1 = "sparse_numeric", e2 = "sparse_numeric"),
          function(e1, e2) sparse_mult(e1, e2))

## -------------------------------
## Plot method
## Plots non-zero positions & values for both vectors; overlaps highlighted
## -------------------------------
setMethod("plot",
          signature(x = "sparse_numeric", y = "sparse_numeric"),
          function(x, y, main = "Sparse vectors (non-zeros)", ...) {
            .check_same_len(x, y)
            ## Build numeric axes limits sans dense coercion
            xmin <- 1
            xmax <- max(1L, x@length)
            ymin <- 0
            ymax <- 0
            if (length(x@value)) ymax <- max(ymax, max(x@value))
            if (length(y@value)) ymax <- max(ymax, max(y@value))
            if (length(x@value)) ymin <- min(ymin, min(x@value))
            if (length(y@value)) ymin <- min(ymin, min(y@value))
            if (ymin == ymax) { ymin <- ymin - 1; ymax <- ymax + 1 }
            
            plot(NA, xlim = c(xmin, xmax), ylim = c(ymin, ymax),
                 xlab = "Position (index)", ylab = "Value",
                 main = main, ...)
            
            if (length(x@pos)) {
              points(x@pos, x@value, pch = 16)
            }
            if (length(y@pos)) {
              points(y@pos, y@value, pch = 1)
            }
            
            ## Overlap highlighting: vertical ticks for indices in both
            ov <- intersect(x@pos, y@pos)
            if (length(ov)) {
              xv <- x@value[match(ov, x@pos)]
              yv <- y@value[match(ov, y@pos)]
              segments(ov, xv, ov, yv, lty = 2)
            }
            legend("topleft",
                   legend = c("x (filled)", "y (open)", "overlap line"),
                   pch = c(16, 1, NA), lty = c(NA, NA, 2), bty = "n")
          })

## -------------------------------
## Extra method (existing generic):
## length() for sparse_numeric
## -------------------------------
setMethod("length", "sparse_numeric", function(x) x@length)

## -------------------------------
## (Optional) Convenience creators
## -------------------------------

## Create from base numeric without as()
sparse <- function(x) as(x, "sparse_numeric")

## Back to base numeric without as()
dense  <- function(x) as(x, "numeric")

