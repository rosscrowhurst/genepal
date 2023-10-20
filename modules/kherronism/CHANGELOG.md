## Source

- Repo: https://github.com/kherronism/rewarewaannotation/tree/1a39a83e22fe2d8665a8c6dc49772cce6579983f
- License: See LICENSE file

## Changes

### repeatmasker

1. Added stub
2. Added author in meta.yml
3. Changed input "tuple val(meta), path(lib)" to "path(lib)"

### braker3

1. Added stub
2. Added author in meta.yml
3. Made output hintsfile optional as it is not produced for ab-initio annotation.
4. Directed `--AUGUSTUS_CONFIG_PATH` to work folder. This avoids "species already exists" error on subsequent runs with same species.
5. Updated version extractor.
6. Added `containerOptions "-B $TMPDIR:$TMPDIR"`