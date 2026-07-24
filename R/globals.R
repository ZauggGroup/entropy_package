# Package-private lookup tables populated lazily by the entropy entry points.
.entropy_cache <- new.env(parent = emptyenv())
