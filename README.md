# PI APROXIMATOR

This program aproximates pi using two old formulas and then compares their results to validate each others correcness.

Compile it with the following command:

```
gcc -g ./pi_record.c -o pi_record -lmpfr -lgmp
```

And then run it passing the precision you want as a command line argument, for example:

```
./pi_record 99999
```

The goal was to beat the 1949 ENIAC (at the time) record of aproximating 2037 decimal places, and by running this program with 1000000 (1 million) as input we can correctly aproximate 301330 decimal places of pi, so mission accomplished.
