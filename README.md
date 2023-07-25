# NL

## installation

1. in SuperCollider run:

```
if(File.exists(Platform.userExtensionDir).not, {File.mkdir(Platform.userExtensionDir)});
```
2. open folder
`Platform.userExtensionDir.openOS;`

3. place the downloaded folder named ``nl`` there

4. recompile class library


## build from source
```shell
# download supercollider source code
cd sc_nl
mkdir build
cd build
# build universal
# change SC_PATH accordingly
cmake -DSC_PATH=/Users/yota/supercollider/ -DSUPERNOVA=ON -DCMAKE_OSX_ARCHITECTURES="x86_64;arm64" ..
cmake -DCMAKE_BUILD_TYPE=RELEASE ..
make
```
