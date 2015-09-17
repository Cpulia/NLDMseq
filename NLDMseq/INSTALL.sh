gcc -fPIC  -g -fPIC  -fPIC   -I/usr/local/include  -c ./SRC/dmatrix.c -o ./SRC/dmatrix.o
gcc -fPIC  -g -fPIC  -fPIC   -I/usr/local/include  -c ./SRC/feature.c -o ./SRC/feature.o
gcc -fPIC  -g -fPIC  -fPIC   -I/usr/local/include  -c ./SRC/gamma.c -o ./SRC/gamma.o
gcc -fPIC  -g -fPIC  -fPIC   -I/usr/local/include  -c ./SRC/gene_expression.c -o ./SRC/gene_expression.o
gcc -fPIC  -g -fPIC  -fPIC   -I/usr/local/include  -c ./SRC/lda.c -o ./SRC/lda.o
gcc -fPIC  -g -fPIC  -fPIC   -I/usr/local/include  -c ./SRC/learn.c -o ./SRC/learn.o        
gcc -fPIC  -g -fPIC  -fPIC   -I/usr/local/include  -c ./SRC/likelihood.c -o ./SRC/likelihood.o
gcc -fPIC  -g -fPIC  -fPIC   -I/usr/local/include  -c ./SRC/newton.c -o ./SRC/newton.o
gcc -fPIC  -g -fPIC  -fPIC   -I/usr/local/include  -c ./SRC/read_beta.c -o ./SRC/read_beta.o
gcc -fPIC  -g -fPIC  -fPIC   -I/usr/local/include  -c ./SRC/util.c -o ./SRC/util.o
gcc -fPIC  -g -fPIC  -fPIC   -I/usr/local/include  -c ./SRC/vbem.c -o ./SRC/vbem.o
gcc -fPIC  -g -fPIC  -fPIC   -I/usr/local/include  -c ./SRC/writer.c -o ./SRC/writer.o
g++ -shared -L/usr/local/lib ./SRC/dmatrix.o ./SRC/feature.o ./SRC/gamma.o ./SRC/gene_expression.o ./SRC/lda.o ./SRC/learn.o ./SRC/likelihood.o ./SRC/newton.o ./SRC/read_beta.o ./SRC/util.o ./SRC/vbem.o ./SRC/writer.o   -o libgene_expre.so  /usr/local/lib/libgsl.a 
