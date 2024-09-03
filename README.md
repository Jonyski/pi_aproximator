# PI APROXIMATOR

This program aproximates pi using two old formulas and then compares their results to validate each others correcness.

Compile it with the following command:

```
gcc -g ./pi_record.c -o pi_record -lmpfr -lgmp
```

And then run it passing the precision you want as a command line argument, for example:

```
./pi_record 9999
```