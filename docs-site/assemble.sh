#!/bin/bash
# Assemble the MkDocs source tree in build/site_src from the repository.
# The tree mirrors the repository layout (README.md, docs/, examples/), so the
# relative links between the pages keep working on the site.
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
SRC="$ROOT/build/site_src"
rm -rf "$SRC"; mkdir -p "$SRC"
cp "$ROOT/docs-site/home.md" "$SRC/index.md"
cp -r "$ROOT/docs-site/assets" "$SRC/assets"
# The repository README becomes the quickstart page; index.md is the gallery.
cp "$ROOT/README.md" "$SRC/quickstart.md"
rsync -a --exclude 'CMakeLists.txt' --exclude 'Doxyfile.in' --exclude 'main.cpp' \
      --exclude 'requirements.txt' "$ROOT/docs/" "$SRC/docs/"
# Examples: pages, figures and the small input files they link to.
rsync -a --exclude '*.TSHS' --exclude '*.psml' --exclude '*.ion.xml' \
      --exclude '__pycache__' "$ROOT/examples/" "$SRC/examples/"
# Links to the repository README point to the quickstart page on the site.
grep -rlE "\]\((\.\./)+README\.md" "$SRC/docs" "$SRC/examples" --include=*.md 2>/dev/null \
  | xargs -r sed -i -E "s#\]\(((\.\./)+)README\.md#](\1quickstart.md#g" || true
echo "assembled $SRC"
