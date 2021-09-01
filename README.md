Simplex
==========

Library implementing simplex algorithm used in operational research

Continuous integration with [Travis-Ci](https://travis-ci.com/quicky2000/simplex) : ![Build Status](https://travis-ci.com/quicky2000/simplex.svg?branch=master)

License
-------
Please see [LICENSE](LICENSE) for info on the license.

Build
-----

Build process is the same used in [Travis file](.travis.yml)
Reference build can be found [here](https://travis-ci.com/quicky2000/simplex)

```
MY_LOCATION=`pwd`
mkdir ../repositories
cd ..
mv $MY_LOCATION repositories
QUICKY_REPOSITORY=`pwd`/repositories
export QUICKY_REPOSITORY
MY_LOCATION=`pwd`
cd $MY_LOCATION/repositories
git clone https://github.com/quicky2000/quicky_tools.git
git clone https://github.com/quicky2000/quicky_exception.git
git clone https://github.com/quicky2000/quicky_utils.git
cd quicky_tools/setup
. setup.sh
cd $MY_LOCATION
chmod a+x repositories/quicky_tools/bin/*
mkdir build
cd build
generate_makefile simplex 
make
```



