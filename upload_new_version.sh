#! /bin/bash

version="$1"

# 
echo "Is everything commited [Y|n]? (git status will be executed)"
git status
read ok

if [[ "$ok" == 'n' || "$ok" == 'N' || "$ok" == 'No' || "$ok" == 'NO' ]]; 
then
    exit
else
    echo "Proceding with $version version upload to PyPI."
fi

ver=$(cat EMDA/_version.py | grep ^__version__ =)
ver_=${ver#*=}

sed -i "" "s|$ver_py_|\'$version\'|g" EMDA/_version.py

#exit

git checkout main
git merge building
git tag $version
git push origin --tags

#gh release create $version
python setup.py sdist
twine upload dist/*



