# Written by Saul Gonzalez Resines
# Mail:sgr40@bath.ac.uk
import warnings
import numpy as np

class DObject:
    def __init__(self,owner=None):
        self.owner=owner
        pass
    # Operator overloading over the class DObject
    # for the sum and the multiplication.

    def _wrapper(self,operator,other):
        # Just create the object
        if self.owner is None:
            print("The variable with name "+self.name+" does not have any owner")
            return None;
        else:
            if isinstance(other,DObject):
                return self.owner.set_ope(operator,self,other);
            elif isinstance(other,float) or isinstance(other,int):
                constant = self.owner.set_cte(other);
                return self.owner.set_ope(operator,self,constant)
            else:
                raise TypeError("Incompatible Types");

    def __add__(self,other):
        return self._wrapper(Add,other);

    def __mul__(self,other):
        return self._wrapper(Multiply,other);

    def __pow__(self,other):
        return self._wrapper(Power,other);
    
    def __div__(self,other):
        return self._wrapper(Powe,other);
    
    def sin(self):
        return self._wrapper(Sine);

class DTrack:
    # This is the object the get all the DObjects
    def __init__(self,head_node=None):
        self.head_node=head_node;
        self.cte={};
        self.var={};
        self.ope={};
        # Numbering the different  operators ad DObjects
        self.n_cte  = 0;
        self.n_var  = 0;
        self.n_add  = 0;
        self.n_mul  = 0;
        self.n_pow  = 0;
        self.n_sine = 0;
        self.n_div  = 0;

    def __repr__(self):
        return "DTrack Object"

    def set_var(self,value,name=None):
        # Set Variable
        this = Variable(self,value,name);
        return this;
        
    def set_cte(self,value,name=None):
        # Set Constant
        this = Constant(self,value,name);
        return this;

       
    def set_ope(self,operator,a,b):
        # Set Operator
        this=operator(self,a,b);
        return this;

    def set_header(self,head):
        # Set the function or header
        self.head_node=head;

    def get_var(self,var):
        if isinstance(Variable):
            return self.var[var.name]
        elif isinstance(str):
            return self.var[var]
        else:
            raise TypeError("Not a Variable or string");
        return None




    def reset(self):
        # Restart the full of the instances (variables)
        self.head_node=None;
        self.cte.clear();
        self.var.clear();
        self.ope.clear();
        self.n_cte=0;
        self.n_var=0;
        self.n_add=0;
        self.n_mul=0;
        self.n_pow=0;
        self.n_div=0;
        self.n_sine=0;

        
    def Topsort(self):
        vis = set();
        order = [];
        def _check(node):
            if node not in vis:
                # avoid the repeat of the same variable or constant;
                # operators are going always to be different;
                vis.add(node);
                if isinstance(node,Operator):
                    for inp in node.inputs:
                        _check(inp);
                order.append(node);

        if self.head_node is None:
            print("The head_node is None so not topolgical order executed");
            return None
        else:
            _check(self.head_node)
            return order;
    
    def Forward(self,order):
        for node in order:
            # Add the placeholder option for an easier handling
            if isinstance(node,Operator):
                node.value = node.forward(*[inp_node.value for inp_node in node.inputs]);

    
    def Backward(self,order):
        # The last of the head_node gradient is always one
        order[-1].gradient=1;
        vis=set();
        for node in reversed(order):
            if isinstance(node,Operator):
                inputs = node.inputs;
                grads  = node.backward(*[inp.value for inp in inputs],dout=node.gradient)
                for inp,grad in zip(inputs,grads):
                    if inp not in vis:
                        inp.gradient = grad;
                    else:
                        inp.gradient += grad;
                vis.add(inp)
        return [node.gradient for node in order]

    def GetGrad(self):
        pass
        
        

class Constant(DObject):
    def __init__(self,owner,value,name=None):
        super().__init__(owner);
        self.owner.n_cte+=1;
        if name is None:
            self.name="cte_"+str(self.owner.n_cte);
        else:
            if name in list(self.owner.cte.keys()):
                self.name="cte_"+str(self.owner.n_cte);
                warnings.warn("Warning the cte with name: "+name+" has been already defined, replacing with "+self.name)
            else:
                self.name=name;
        self.value=value;
        # Updating the owner dictionary
        # and assign the object to the dictionary position
        self.owner.cte[self.name]=self;
        self=self.owner.cte[self.name];


    def __repr__(self):
        return "The cte is : name "+self.name+" with value "+str(self.value)

class Variable(DObject):
    def __init__(self,owner,value,name=None):
        super().__init__(owner);
        self.owner.n_var+=1;
        if name is None:
            self.name="var_"+str(self.owner.n_var);
        else:
            if name in list(self.owner.var.keys()):
                self.name="cte_"+str(self.owner.n_var);
                warnings.warn("Warning the variable with name: "+name+" has been already defined, replacing with "+self.name)
            else:
                self.name=name;
        self.value=value;
        # Updating the owner dictionary
        # and assign the object to the dictionary position
        self.owner.var[self.name]=self;
        self=self.owner.var[self.name];

    def __repr__(self):
        return "The variable is: name "+self.name+" with value "+str(self.value);


class Operator(DObject):
    def __init__(self,owner,name="ope_"):
        super().__init__(owner);
        self.value=None;
        self.inputs=[];
        self.name=name;
        self.gradient=None;


    def __repr__(self):
        return "The name of the operator is, "+self.name ;

class Add(Operator):
    def __init__(self,owner,a,b,name="add_"):
        super().__init__(owner,name);
        self.owner.n_add+=1;
        self.name="add_"+str(self.owner.n_add)
        self.inputs=[a,b]
        # Add the operator to the owner
        self.owner.ope[self.name]=self;
        self=self.owner.ope[self.name]

    def forward(self,a,b):
        return a+b;

    def backward(self,a,b,dout):
        return dout,dout;

class Multiply(Operator):
    def __init__(self,owner,a,b,name="mul_"):
        super().__init__(owner,name);
        self.owner.n_mul+=1;
        self.name="mul_"+str(self.owner.n_mul);
        self.inputs=[a,b]
        # Add the operator to the owner
        self.owner.ope[self.name]=self;
        self=self.owner.ope[self.name];


    def forward(self,a,b):
        return a*b;

    def backward(self,a,b,dout):
        return dout*b,dout*a;

class Power(Operator):
    def __init__(self,owner,a,b,name="pow_"):
        super().__init__(owner,name);
        self.owner.n_pow+=1;
        self.name="pow_"+str(self.owner.n_mul);
        self.inputs=[a,b]
        # Add the operator to the owner
        self.owner.ope[self.name]=self;
        self=self.owner.ope[self.name];

    def forward(self,a,b):
        return np.power(a,b)

    def backward(self,a,b,dout):
        return dout*b*np.power(a,b-1.0),dout*np.log(a)*np.power(a,b)

class Divide(Operator):
    def __init__(self,owner,a,b,name="div_"):
        super().__init__(owner,name);
        self.owner.n_div+=1;
        self.name="div_"+str(self.owner.n_div);
        self.inputs=[a,b]
        # Add the operator to the owner
        self.owner.ope[self.name]=self;
        self=self.owner.ope[self.name];

    def forward(self,a,b):
        return a/b;

    #def backward(self,a,b):
    #    return dout*(),dout*

class Sine(Operator):
    def __init__(self,owner,a,name="sine_"):
        super().__init__(owner,name);
        self.owner.n_sine+=1;
        self.name="sine_"+str(self.owner.n_sine);
        self.inputs=[a]
        # Add the operator to the owner
        self.owner.ope[self.name]=self;
        self=self.owner.ope[self.name];
    
    def forward(self,a):
        return np.sin(a);

    def backward(self,a,dout):
        return dout*np.cos(a);


class Cosine(Operator):
    # Here the Cosine
    pass

class Tangent(Operator):
    # Here the Tangent
    pass





if __name__=='__main__':
    import sys,os;
    '''
    Just make the autodifferentiation
    with the header that you are looking.
    '''
    myform=DTrack();
    x  = myform.set_var(value=5.0,name="x");
    c1 = myform.set_cte(value=2.0);
    c2 = myform.set_cte(value=3.0);
    header = x**2
    # Maybe the order should be an internal variable
    myform.set_header(header);
    order = myform.Topsort();
    myform.Forward(order); 
    grads=myform.Backward(order);
    myform.reset()
    
    
    
