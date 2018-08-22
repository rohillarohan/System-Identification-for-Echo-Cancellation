# Adaptive Algorithms for Echo Cancellation
    
System Identification/Modeling for Echo Cancellation:   In this application an adaptive filter is used to identify the impulse response of the path between the source from which the echo originates adn the point where it appears. The output of the adaptive filter, which is an estimate of the echo signal, can be used to cancel the undesirable echo. 

In my Project I simulated a modeling problem by feeding unit variance white Gaussian noise through a system H(z) (which can be considered a as a channel) the output from this system is fed into the plant we want to estimated G(z). The plant output is then contaminated with an additive white Gaussian Noise with variance sigma^2. The output after the AWGN is considered a desired signal d(n). An N-tap adaptive filter W(z), is used to estimate the plant G(z). 

We estimate the plant G(z) by minimizing the mean square error (MSE) between the output from the adaptive filter W(z) which is y(n) and d(n).

To minimize this cost function we have used different techniques:
    1) Pure LMS Algorithm
    2) Sign LMS Algorithm
    3) Sign-Regressor LMS Algorithm
    4) Sign-Sign LMS Algorithm
    5) Normalized LMS Algorithm
    6) Affine Projection LMS Algorithm
    7) Recursive Least Square (RLS) Algorithm

At the end all the algorithms are compared.
