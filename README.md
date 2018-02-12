### Deterministic Identification of Discrete-Time State Space Systems.

[Ai,Bi,Ci,Di] = deti(u,y)  

Function deti(u,y) identifies a discrete-time state space system (calculates matrices A,B,C,D and system's order) using input/output measurements only. The function does not require large amount of input/output measurements in order to produce good and stable results. An amount of about 45-50 measurements is usually enough.

Input argument u: Input measurements of the system we want to identify.
Input argument y: Output measurements of the system we want to identify.

The function returns the identified matrices Ai,Bi,Ci,Di, the rank of the system, and the goodness of fit between the identified and the original system.

There is a small chance for the function to produce a system with deficient order(rank). Depending on the amount of input/output measurements the function will produce a better estimation.

Note that this function only works for (linear) discrete-time state space systems, without disturbances in the input/output measurements.

### References

[1] Moonen M., DeMoor B., Vandenberghe L., Vandewalle J., On And Off-line
Identification of Linear State Space Models, International Journal of Control,
Vol. 49(1), pp. 219-232, 1989
