import numpy as np
import matplotlib.pyplot as plt
import sys



def fun(x,y,l):
    return ((x-2.0)**2)+((y-5.0)**2)+l*(x+y-8.0) ;

def gradx(x,y,l):
    return 2*(x-2.0)+l

def grady(x,y,l):
    return 2*(y-5.0)+l

def gradl(x,y,z):
    return -1*(x+y-8.0)




def GradientDescent(start,step,fun,max_steps,tol=0.0001,beta=0.8,grads=[gradx,grady,gradl]):
    # Gradient Descent computing the derivative with finite differences
    # The h of the derivative is a really small number to get the error to the minimum.
    hdif=0.001;
    # Starting tolerance with a big number so that we will not miss the 
    st_tol=1000.0;
    print("Maxsteps,",max_steps)
    # Start with = iterations
    count=0;
    # Check step
    if type(step) == int or float:
        tmp=step;
        step=np.zeros(len(start));
        step[:]=tmp;
    elif type(step)==list:
        step=np.array(step);
    else:
        pass

    if len(step)==len(start):
        pass
    else:
        print("Step does not have the same N of parameters as start");
        return ;
    or_step=step.copy();
    x=np.array(start)
    xbef=x;
    while (st_tol>=tol):
        step=or_step.copy();
        der=[];
        # Evaluate the function at the actual iteration
        actual=fun(*x)
        for c,i in enumerate(x):
            # Derivative computation (Finite differences)
            x1=x.copy()
            x1[c]=x[c]+hdif;
            derivative = (fun(*x1)-actual)/hdif
            derivative1 = grads[c](*x)
            print("compare derivatives,",derivative,derivative1)
            der.append(derivative1)
            ## Backtracking line search
            #x1[c] = x[c]-(step[c]*derivative);
            #its=0;
            ## Back Tracking Line Search
            #while (fun(*x1) > ( actual - ((step[c]/2.0)*abs(derivative)**2)) ):
            #    # Update the step 
            #    step[c]=step[c]*beta;
            #    # Update the x1
            #    x1[c]=x[c]-step[c]*derivative
            #    its=its+1;
        print("The step is,",step)
        print(der)
        # Here update the X, the count and the step
        der=np.array(der);
        xbef=x.copy();
        x=x-(step*der);
        count=count+1;
        #st_tol=np.linalg.norm(xbef-x)
        print(st_tol)
        if count>=max_steps:
            print("Maximum steps reached")
            break;
    return x
            

step=0.01;
start=[6.0,7.0,-2.0]
max_steps=3000
tol=0.0001
where=GradientDescent(start,step,fun,max_steps,tol=tol);
print("result,",where)
print("sum",where[0]+where[1])

