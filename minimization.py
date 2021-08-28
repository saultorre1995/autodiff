# Written by Saul Gonzalez Resines
# Mail:sgr40@bath.ac.uk
# Bibliography : Algorithms for Optimization Mykel J. Kochenderfer Tim A. Wheeler
import autodiff as au
import numpy as np
import sys
import tracemalloc
# Track memory


class Minimization:
    
    '''
    Minimization algorithm 
    for the graphs created in the 
    autodifferentiation package
    '''

    def __init__(self,dtrack,variables=None,constraints=None):
        # Graph for the Minimization (with all the variables for the minimization)
        self.dtrack = dtrack;
        if variables is None:
            pass
        # Set the starting value for the variable
        elif type(variables)==dict:
            for varname in list(varibles.keys()):
                self.dtraj.var[varname].value=variables[varname];
        # Set the constraints for the calculation of the minimum
        self.constraints = constraints;
        self.output={};

class Gradescent(Minimization):
    '''
    Gradescent: Unconstrained algorithm
    '''
    def __init__(self,dtrack,variables=None,constraints=None):
        super().__init__(dtrack,variables,constraints);

    def run(self,steps,alpha=0.1):
        varnames = self.dtrack.var.keys();
        order = self.dtrack.Topsort();
        mytime=0
        for i in range(steps):
            # Calculate the gradient
            self.dtrack.Forward(order);
            self.dtrack.Backward(order);
            # Calculate the gradient of the var
            gradients = np.array([self.dtrack.var[varname].gradient for varname in varnames])
            # Normalize and scale the gradients
            gradients = (gradients/np.linalg.norm(gradients))*alpha
            # Assign the new values for the variables
            for c,varname in enumerate(varnames):
                self.dtrack.var[varname].value=self.dtrack.var[varname].value -gradients[c]
        # The return value is dictionary of the values of the minimums of the variables
        self.output = {variable.name:variable.value  for variable in self.dtrack.var.values()}

class Adagrad(Minimization):
    '''
    Adagrad minimization algorithm
    Unconstrained algorithm
    '''
    def run(self,steps,alpha=0.01):
        # eps avoid the multiplication by zero
        # alpha is the learning rate, normally 0.01
        eps = 10e-10;
        varnames = self.dtrack.var.keys();
        order = self.dtrack.Topsort();
        # Accumulative vector
        vector = np.zeros(len(varnames))
        for i in range(steps):
            self.dtrack.Forward(order);
            self.dtrack.Backward(order);
            # Calculate the gradient of the var
            gradients = np.array([self.dtrack.var[varname].gradient for varname in varnames])
            vector = vector + np.power(gradients,2);
            # Update values
            for c,varname in enumerate(varnames):
                self.dtrack.var[varname].value=self.dtrack.var[varname].value - gradients[c]*( alpha/(eps+np.sqrt(vector[c])))
        self.output = {variable.name:variable.value  for variable in self.dtrack.var.values()}


class Nesterov(Minimization):
    '''
    Nesterov minimization momentum
    Unconstrained algorithm
    '''
    def run(self,steps,alpha=0.01,beta=0.1):
        varnames = self.dtrack.var.keys();
        order = self.dtrack.Topsort();
        # start with momentum equal to zero
        mom = 0
        for i in  range(steps):
            self.dtrack.Forward(order);
            self.dtrack.Backward(order);
            gradients = np.array([self.dtrack.var[varname].gradient for varname in varnames])
            # Calculate momentum
            mom = beta * mom - alpha * gradients
            for c,varname in enumerate(varnames):
                self.dtrack.var[varname].value=self.dtrack.var[varname].value + mom[c]
        # get the output
        self.output = {variable.name:variable.value  for variable in self.dtrack.var.values()}



class Adadelta(Minimization):
    '''
    Adadelta Unconstrained Minimization
    '''
    def run(self,steps,beta_grad=0.9,beta_diff=0.9,eps=10e-8):
        varnames = self.dtrack.var.keys();
        order = self.dtrack.Topsort();
        # Initialize the gradient squared sum and the difference squared sums
        grad_sq_sum  = np.zeros(len(varnames));
        diff_sq_sum = np.zeros(len(varnames));
        for i in  range(steps):
            self.dtrack.Forward(order);
            self.dtrack.Backward(order);
            gradients = np.array([self.dtrack.var[varname].gradient for varname in varnames])
            # update grad_sq_sum
            grad_sq_sum = beta_grad*grad_sq_sum + (1.0-beta_grad)*gradients*gradients;
            # calculate delta_x that start in 0.0
            delta_x = - ((np.sqrt(diff_sq_sum) + eps)/(eps+np.sqrt(grad_sq_sum)))*gradients
            #print(delta_x)
            diff_sq_sum = beta_diff*diff_sq_sum + (1.0 - beta_diff) * delta_x * delta_x
            for c,varname in enumerate(varnames):
                self.dtrack.var[varname].value=self.dtrack.var[varname].value + delta_x[c]
        # get the output
        self.output = {variable.name:variable.value  for variable in self.dtrack.var.values()}


class Adam(Minimization):
    '''
    Adam Unconstrained Minimization
    '''
    def run(self):
        
        pass
 


