#include <stdio.h>
#include <math.h>

/* generate a table of values for the funtion

   10 * lg(1 + 2^(x/10))

   where lg() is the binary log. */

int
main () {
    int    j;
    double i, result, scale, delta, max, lg2;

    scale = 100;
    delta = 1 / scale;
    max = -100.0;
    lg2 = log(2);

    printf("float zFloatAddScale = %15.32f;\n", scale);
    printf("float zFloatAddMax = %15.32f;\n", max);
    printf("float zFloatAddTable[]={\n");

    j = 0;
    i = -j / scale;
    result = i * lg2;
    result = 10*log(1+exp(result / 10 * lg2));
    result /= lg2;
    printf ("%15.32f", (float)result);
    j++;

    while (i > max) {
	result = i;
	result = 10*log(1+exp(result / 10 * lg2));
	result /= lg2;
	printf (",\n%15.32f", (float)result);

	j++;
	i = -j / scale;
    }
    printf("};\n");

    printf ("\nint zFloatAddCount = %d;\n", j);

    return 0;
}
