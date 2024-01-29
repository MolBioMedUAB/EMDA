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

ver_py=$(cat setup.py | grep ^version=)
ver_py_=${ver_py#*=}

ver_cfg=$(cat setup.cfg | grep '^version')
ver_cfg_=${ver_cfg#*=}

echo $ver_py_
echo $ver_cfg_

sed -i "" "s|$ver_py_|\'$version\'|g" setup.py
sed -i "" "s|$ver_cfg_| $version|g" setup.cfg

#exit

#git checkout main
#git merge building
git tag $version
git push origin --tags

#gh release create $version
python setup.py sdist
twine upload dist/*



