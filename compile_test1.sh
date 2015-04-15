g++ test1.cc -o test1.exe `root-config --cflags --libs` -I$LCIO/include -L$LCIO/lib -llcio -llcioDict -I$LCIO_READ/include -L$LCIO_READ/lib -llcioread
#g++ test1.cc -o test1.exe `root-config --cflags --libs` -I$LCIO/include -L$LCIO/lib -llcio -llcioDict -I/home/rafopar/WORKDIR/HPS/LCIO_Read/include -L/home/rafopar/WORKDIR/HPS/LCIO_Read/lib -llcioread
