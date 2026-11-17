#!/bin/bash
# prune_keep.sh ROOT KEEPLIST
#   keep only the keep-list's entries of the package tree ROOT; the rest is
#   deleted. Keep-list lines, relative to ROOT ('#' starts a comment):
#     d <dir>    the directory's files and symlinks (not its subdirectories)
#     f <file>   one file or symlink
#   A listed entry missing from ROOT stops the build: the list no longer
#   matches the package.
set -euo pipefail
root=$1
list=$2
new=${root}.keep
rm -rf "$new"
mkdir -p "$new"
cd "$root"
grep -v '^[[:space:]]*\(#\|$\)' "$list" | while read -r type path; do
  case $type in
    d) test -d "$path"
       find "$path" -mindepth 1 -maxdepth 1 ! -type d -print0 \
         | xargs -0 -r cp -a --parents -t "$new" ;;
    f) cp -a --parents "$path" "$new" ;;
    *) echo "prune_keep.sh: bad line: $type $path" >&2; exit 1 ;;
  esac
done
cd /
before=$(du -sm "$root" | cut -f1)
rm -rf "$root"
mv "$new" "$root"
echo "prune_keep.sh: $root ${before} MB -> $(du -sm "$root" | cut -f1) MB"
