#! /bin/bash

version="$1"

if [[ "$version" == '' ]];
then
    echo "Please, add the version as an argument"
    exit
fi

# 
echo "Is everything commited [Y|n]? (git status will be executed)"
git status
read ok

if [[ "$ok" == 'y' || "$ok" == 'Y' || "$ok" == 'Yes' || "$ok" == 'YES' || "$ok" == '' ]]; 
then
    echo "Proceding with $version version upload to PyPI."
else
    exit
fi

ver=$(cat EMDA/_version.py | grep '^__version__ =')
ver_=${ver#*=}
echo $ver_

sed -i "" "s|$ver_|\'$version\'|g" EMDA/_version.py

#exit

git commit EMDA/_version.py -m "Update to version $version"

black EMDA/
branch=$(git rev-parse --abbrev-ref HEAD)

if [[ "$branch" != 'main' ]];
then
    git checkout main
    git merge building
fi

git tag $version
git push origin --tags

#gh release create $version
python setup.py sdist
twine upload dist/*



