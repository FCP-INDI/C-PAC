#!/usr/bin/env bash

# Copyright (C) 2024-2025  C-PAC Developers

# This file is part of C-PAC.

# C-PAC is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.

# C-PAC is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
# License for more details.

# You should have received a copy of the GNU Lesser General Public
# License along with C-PAC. If not, see <https://www.gnu.org/licenses/>.


set -euo pipefail
trap 'echo "❌ Script failed at line $LINENO with exit code $?"' ERR

# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------

git_add_with_retry() {
  local file=$1
  local attempts=0
  local max_attempts=10
  while ! git add "$file"; do
    attempts=$((attempts+1))
    echo "Git add failed for $file (attempt $attempts), retrying..."
    sleep 1
    if [[ $attempts -ge $max_attempts ]]; then
      echo "❌ Failed to git add $file after $max_attempts attempts"
      exit 1
    fi
  done
}

update_file_if_changed() {
  # Run a regex replacement or copy on a file and stage it if it changed
  local expr=$1
  local src=$2
  local dest=${3:-$src}

  local changed=0
  if [[ -n "$expr" ]]; then
    tmp=$(mktemp)
    sed -E "$expr" "$src" > "$tmp"
    if ! cmp -s "$tmp" "$dest"; then
      mv "$tmp" "$dest"
      git_add_with_retry "$dest"
      changed=1
    else
      rm "$tmp"
    fi
  else
    if [[ ! -f "$dest" ]] || ! cmp -s "$src" "$dest"; then
      cp "$src" "$dest"
      git_add_with_retry "$dest"
      changed=1
    fi
  fi
  return $changed
}

log_info() {
  echo "=== $* ==="
}

# -------------------------------------------------------------------------
# Main
# -------------------------------------------------------------------------

START_DIR=$(pwd)
SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
REPO_ROOT="$(realpath "$SCRIPT_DIR/../..")"

# -------------------------------------------------------------------------
# Fetch version
# -------------------------------------------------------------------------
log_info "Fetching version"
VERSION=$(python -c "import sys; sys.path.insert(0, '$REPO_ROOT/CPAC'); from info import __version__; print(__version__.split('+', 1)[0])")
VERSION_FILE="$REPO_ROOT/version"
if [[ -f "$VERSION_FILE" ]]; then
    cd "$REPO_ROOT"
    OLD_VERSION=$(git show "$(git log --pretty=format:'%h' -n 1 -- version | tail -n 1)":version)
    cd "$START_DIR"
else
    OLD_VERSION="<none>"
fi
echo "v${VERSION}" > "$VERSION_FILE"

# -------------------------------------------------------------------------
# Write version file and stage it
# -------------------------------------------------------------------------
log_info "Updating version file"
if update_file_if_changed "" <(echo "v${VERSION}") "$VERSION_FILE"; then
  git_add_with_retry "$VERSION_FILE"
fi

# -------------------------------------------------------------------------
# Update YAML config files
# -------------------------------------------------------------------------
log_info "Updating YAML config files"
VERSION_EXPR="s/^(# [Vv]ersion ).*$/# Version ${VERSION}/g"
for YAML_FILE in "$REPO_ROOT"/CPAC/resources/configs/{*.yml,*.yaml,test_configs/*.yml,test_configs/*.yaml}; do
  [[ -e "$YAML_FILE" ]] || continue

  echo "Processing ${YAML_FILE}"
  echo "Applying regex: ${VERSION_EXPR}"

  # Run sed safely
  tmp=$(mktemp)
  if ! sed -E "$VERSION_EXPR" "$YAML_FILE" > "$tmp"; then
    echo "❌ sed failed on $YAML_FILE"
    rm "$tmp"
    exit 1
  fi

  if ! cmp -s "$tmp" "$YAML_FILE"; then
    mv "$tmp" "$YAML_FILE"
    echo "Updated $YAML_FILE"
    git_add_with_retry "$YAML_FILE"
  else
    rm "$tmp"
    echo "No changes needed for $YAML_FILE"
  fi
done

# -------------------------------------------------------------------------
# Update Dockerfiles (only C-PAC tags)
# -------------------------------------------------------------------------
log_info "Updating Dockerfiles"
NEW_VERSION=$(<"$VERSION_FILE")

if [[ "$OLD_VERSION" != "$NEW_VERSION" ]]; then
  for DOCKERFILE in "$REPO_ROOT"/.github/Dockerfiles/*.Dockerfile; do
    if grep -q "FROM ghcr\.io/fcp-indi/c-pac/.*-${OLD_VERSION}" "$DOCKERFILE"; then
      echo "Updating C-PAC version in ${DOCKERFILE} from ${OLD_VERSION} to ${NEW_VERSION}"

      if [[ "$OSTYPE" == "darwin"* ]]; then
        # macOS sed
        sed -i "" "s/-${OLD_VERSION}/-${NEW_VERSION}/g" "$DOCKERFILE"
      else
        # Linux sed
        sed -i -E "s/-${OLD_VERSION}/-${NEW_VERSION}/g" "$DOCKERFILE"
      fi

      git_add_with_retry "$DOCKERFILE"
    fi
  done
fi

# -------------------------------------------------------------------------
# Overwrite top-level Dockerfiles
# -------------------------------------------------------------------------
log_info "Updating top-level Dockerfiles"
TOP_DOCKERFILES=(
  ".github/Dockerfiles/C-PAC.develop-jammy.Dockerfile:Dockerfile"
  ".github/Dockerfiles/C-PAC.develop-lite-jammy.Dockerfile:variant-lite.Dockerfile"
)
for SRC_DST in "${TOP_DOCKERFILES[@]}"; do
  # Split SRC_DST by colon safely
  SRC="${SRC_DST%%:*}"
  DST="${SRC_DST##*:}"

  FULL_SRC="$REPO_ROOT/$SRC"
  FULL_DST="$REPO_ROOT/$DST"

  if [[ ! -f "$FULL_SRC" ]]; then
    echo "⚠️ Source Dockerfile does not exist: $FULL_SRC"
    continue
  fi
  echo "Updating top-level Dockerfile: $FULL_DST from $FULL_SRC"
  cp "$FULL_SRC" "$FULL_DST" && git_add_with_retry "$FULL_DST"
done

# Return to original directory
cd "$START_DIR"

# -------------------------------------------------------------------------
# Summary
# -------------------------------------------------------------------------
echo
echo "Version changed: (from ${OLD_VERSION} to ${NEW_VERSION})"
echo "======================"
