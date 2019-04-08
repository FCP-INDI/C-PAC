#!/bin/bash

VERSION=$1

VERSION_PIECES=(${VERSION//./ })

if [[ ! "$VERSION" =~ ^[0-9]+\.[0-9]+\.[0-9]+(\.[a-zA-Z0-9]+)?$ ]]; then
    echo "Please, provide an argument with the new version in the following format: [0-9].[0-9].[0-9]"
    exit
fi

MAJOR=${VERSION_PIECES[0]}
MINOR=${VERSION_PIECES[1]}
MACRO=${VERSION_PIECES[2]}
EXTRA=${VERSION_PIECES[3]}

echo -n "v$VERSION" > version

sed -Ei "s/^_version_major\s*=\s*[0-9]+$/_version_major = ${MAJOR}/g" CPAC/info.py
sed -Ei "s/^_version_minor\s*=\s*[0-9]+$/_version_minor = ${MINOR}/g" CPAC/info.py
sed -Ei "s/^_version_micro\s*=\s*[0-9]+$/_version_micro = ${MACRO}/g" CPAC/info.py
sed -Ei "s/^_version_extra\s*=\s*[\"'][a-zA-Z0-9]*[\"']$/_version_extra = '${EXTRA}'/g" CPAC/info.py

rename_patterns () {
    sed -Ei "1,5s/^\# Version [0-9]+\.[0-9]+\.[0-9]+(\.[a-zA-Z0-9]+)?$/# Version $1/" $2
}

export -f rename_patterns
find -regex ".*\.ya?ml$" -exec bash -c "rename_patterns "${VERSION}" \"\$0\"" {} \;