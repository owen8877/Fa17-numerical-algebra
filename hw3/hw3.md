## Numerical Linear Algebra - hw3

_陈子恒 1500010632_

### Prob 1

The result shows as follows:

| normal gauss | improved gauss | normal cholesky | improved cholesky |     QR     | QR(modified) |
| :----------: | :------------: | :-------------: | :---------------: | :--------: | :----------: |
|   7.2593e8   |   8.5973e-08   |        -        |         -         |    NaN     |    6.7612    |
|  1.2856e-17  |   1.2856e-17   |   4.2826e-18    |    1.2856e-17     | 2.4289e-17 |  2.4876e-17  |
|   4.7865e2   |    6.4146e2    |    7.4571e7     |     3.1407e2      |  4.0672e2  |   6.0747e2   |

It seems that the householder algorithm provided by the textbook has poor numerical stability.

_Note : The modified QR method uses the following householder algorithm :_

```matlab
% b = 2*v(1)^2/(sigma+v(1)^2);
% v = v / v(1);
b = 1;
v = v / sqrt((sigma+v(1)^2)/2);
```

### Prob 2

```plain
LS Result:
	a	1.000000e+00
	b	1.000000e+00
	c	1.000000e+00
Error	6.080942e-16
```

### Prob 3

```plain
LS Result:
    2.0775
    0.7189
    9.6802
    0.1535
   13.6796
    1.9868
   -0.9582
   -0.4840
   -0.0736
    1.0187
    1.4435
    2.9028

Error	1.985495e-12
```

