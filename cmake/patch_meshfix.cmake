# Idempotently guard MeshFix's closeLoadingSession() against fclose(NULL).
#
# meshfix_from_eigen_matrices() (examples/meshfix_eigen.h) builds the mesh in
# memory and calls closeLoadingSession(NULL, ...), whose final `fclose(fp)`
# segfaults on modern glibc when fp is NULL.  Guarding fclose is a no-op for the
# genuine file-loading callers (where fp is valid).
#
# Invoked as a FetchContent PATCH_COMMAND with the working directory set to the
# MeshFix source tree, so MESHFIX_IO_CPP is a path relative to that tree.
file(READ "${MESHFIX_IO_CPP}" _contents)
# normalize (so re-running the patch is a no-op) then guard
string(REPLACE "if(fp) fclose(fp);" "fclose(fp);" _contents "${_contents}")
string(REPLACE "fclose(fp);" "if(fp) fclose(fp);" _contents "${_contents}")
file(WRITE "${MESHFIX_IO_CPP}" "${_contents}")
