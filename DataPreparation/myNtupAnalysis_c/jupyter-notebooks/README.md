Make sure C++ library is compiled:

```
g++ -shared -fPIC -o Cfunctions.so Cfunctions.cxx `root-config --cflags --glibs`
```