# ODE

TODO
1. how to define init value for gene data, variance?
2. how to define the max iter, n and R in generate trajectories?
2. A: iter for the GPLVM data
3. need explain the last module in cost function(Guo's: suppose different states would appear in the same possibility. 
The last part need to increase or remove).
4. hill coefficient need special definition
5. need some trick to deal inf or nan(add random number for example).
6. The bad combination of params would lead to unstable state more easily, it could be a standard for end up the loop in 
trajectories operation
7. Why choose some genes to GPLVM model, instead of all genes?(too large?)
8. the final result would still have unstable state when test on different gene init value.
9. Attr_param defines method?