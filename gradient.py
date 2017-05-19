from math import *


epslon = pow(10,-10)
ro = 0.01
my_values = ()
iteracoes_maximas = pow(10,4)
margem = pow(10,-10)
def multiplica_matrix_22_22(matrix_a, matrix_b):
    nova_matrix = []
    first_line = [(matrix_a[0][0]*matrix_b[0][0])+(matrix_a[0][1]*matrix_b[1][0]), (matrix_a[0][0]*matrix_b[0][1])+(matrix_a[0][1]*matrix_b[1][1])]
    second_line = [(matrix_a[1][0]*matrix_b[0][0]+matrix_a[1][1]*matrix_b[1][0]), (matrix_a[1][0]*matrix_b[0][1]+matrix_a[1][1]*matrix_b[1][1])]
    nova_matrix.append(first_line)
    nova_matrix.append(second_line)
    return nova_matrix

def multiplica_matrix_12_22(matrix_a,matrix_b):
    return [[(matrix_a[0][0]*matrix_b[0][0] + matrix_a[0][1]*matrix_b[1][0]), matrix_a[0][0]*matrix_b[0][1]+matrix_a[0][1]*matrix_b[1][1]]]

def multiplica_matrix_22_21(matrix_a,matrix_b):
    return [[(matrix_a[0][0] * matrix_b[0]) + (matrix_a[0][1]*matrix_b[1])],[(matrix_a[1][0]*matrix_b[0]) + (matrix_a[1][1]*matrix_b[1])]]

def transposta_12(my_matrix):
    nova_matrix = [[my_matrix[0][0]],[my_matrix[0][1]]]
    return nova_matrix


def inverte_quadrada(matrix):
    det = (matrix[0][0]*matrix[1][1])-(matrix[0][1]*matrix[1][0])
    # det= det * (-1) para se adequar ao metodo
    det = det * (-1)
    matrix_invertida = [[(matrix[1][1]/(1.0*det)),(-matrix[0][1]/(1.0*det))],[(-matrix[1][0]/(1.0*det)),matrix[0][0]/(1.0*det)]]
    return matrix_invertida

def gradiente(x,y,function):
    if (function == 1):
        x_grad = (1-2*x)/((x-1)*x)
        y_grad = (1-2*y)/((y-1)*y)
        return (x_grad,y_grad)
    elif(function == 2):
        x_grad = ((1-2*x)/(2*(x-1)*x*sqrt(-log((x-1)*x*(y-1)*y))))
        y_grad = ((1-2*y)/(2*(y-1)*y*sqrt(-log((x-1)*x*(y-1)*y))))
        return(x_grad,y_grad)
    elif(function== 3):
        x_grad = 2 * x
        y_grad = 2 * y
        return(x_grad,y_grad)

def gradiente_2(function, x,y):
    if (function == 1):
        hessiana = []
        grad_1xx = (2*pow(x,2)-2*x+1)/(pow((x-1),2)*pow(x,2))
        grad_1yy = 0
        grad_2xx = 0
        grad_2yy = (2*pow(y,2)-2*y+1)/(pow((y-1),2)*pow(y,2))
        hessiana.append([grad_1xx,grad_1yy])
        hessiana.append([grad_2xx,grad_2yy])
        return hessiana
    elif(function == 2):
        hessiana = []
        grad_1xx = - ((4*pow(x,2)-4*x+2)*log((x-1)*x*(y-1)*y)+pow((1-2*x),2))/(4*pow((x-1),2)*pow(x,2)*pow((-log((x-1)*x*(y-1)*y)),1.5))
        grad_1yy = - ((2*x-1)*(2*y-1))/(4*(x-1)*x*(y-1)*y*pow((-log((x-1)*x*(y-1)*y)),1.5))
        grad_2xx = - ((2*x-1)*(2*y-1))/(4*(x-1)*x*(y-1)*y*pow((-log((x-1)*x*(y-1)*y)),1.5))
        grad_2yy = - ((4*pow(y,2)-4*y+2)*log((x-1)*x*(y-1)*y)+pow((1-2*y),2))/(4*pow((y-1),2)*pow(y,2)*pow((-log((x-1)*x*(y-1)*y)),1.5))
        hessiana.append([grad_1xx,grad_1yy])
        hessiana.append([grad_2xx,grad_2yy])
        return hessiana

def my_function(function,x,y):
    if (function == 1): return -log(x*(1-x)*y*(1-y))
    if (function == 2): return sqrt(-log(x*(1-x)*y*(1-y)))
    if (function == 3): return pow(x,2) + pow(y,2)

def phi(function,t,d):
    x = my_values[0] + d[0] * t
    y = my_values[1] + d[1] * t
    if (function == 1): return -log(x*(1-x)*y*(1-y))
    if (function == 2): return sqrt(-log(x*(1-x)*y*(1-y)))
    #f_3 e uma funcao simples para teste
    #print("Valor de X: ",x)
    #print("Valor de x a quarta:",pow(x,4))
    #print("Valor de Y: ",y)
    #print("Valor de y a quarta:",pow(y,4))
    if (function == 3): return pow(x,2) + pow(y,2)


def secao_aurea(funcao,epslon,ro,d):
    theta_1 = (3-sqrt(5))/2
    theta_2 = 1 - theta_1
    # obtencao do intervalo [a, b]
    a, s, b = 0, ro, 2*ro
    while (phi(funcao,b,d) < phi(funcao,s,d)):
        a, s, b = s, b, 2*b
    #  obtencao de t*
    u = a + theta_1*(b-a)
    v = a + theta_2*(b-a)
    while((b-a)>epslon):
        if(phi(funcao,u,d) < phi(funcao,v,d)):
            b,v = v, u
            u = a + theta_1*(b-a)
        else:
            a,u= u,v
            v = a + theta_2*(b-a)
    t = (u+v)/2
    return t

def metodo_gradiente(funcao, x_inicial, y_inicial,iteracao = 0):
    global my_values
    my_values = (x_inicial,y_inicial)
    old_values = (0,0)
    while( ( (fabs(old_values[0]-my_values[0]) > margem) and (fabs(old_values[1]-my_values[1]) > margem) ) and  iteracao < iteracoes_maximas ):
        old_values = my_values
        d = (-1 * gradiente(my_values[0],my_values[1],funcao)[0],-1*gradiente(my_values[0],my_values[1],funcao)[1])
        t = secao_aurea(funcao, epslon, ro, d)
        my_values = (my_values[0] + d[0]*t, my_values[1] + d[1]*t)
        iteracao += 1
    return (my_values,iteracao,fabs(old_values[0]-my_values[0]),my_function(funcao,my_values[0],my_values[1]))

def metodo_newton(funcao, x_inicial, y_inicial,iteracao = 0):
    global my_values
    my_values = (x_inicial,y_inicial)
    old_values = (0,0)
    while( ( (fabs(old_values[0]-my_values[0]) > margem) and (fabs(old_values[1]-my_values[1]) > margem) ) and  iteracao < iteracoes_maximas ):
        old_values = my_values
        hessiana = gradiente_2(funcao, my_values[0], my_values[1])
        hessiana_i =  inverte_quadrada(hessiana)
        matrix_d = multiplica_matrix_22_21(hessiana_i,gradiente(my_values[0],my_values[1],funcao))
        d = (matrix_d[0][0],matrix_d[1][0])
        t = secao_aurea(funcao, epslon, ro, d)
        my_values = (my_values[0] + d[0]*t, my_values[1] + d[1]*t)
        iteracao += 1
    return (my_values,iteracao,fabs(old_values[0]-my_values[0]),my_function(funcao,my_values[0],my_values[1]))

def metodo_quase_newton(funcao, x_inicial, y_inicial,iteracao = 0):
    global my_values
    my_values = (x_inicial,y_inicial)
    old_values2 = (0,0)
    while( ( (fabs(old_values2[0]-my_values[0]) > margem) and (fabs(old_values2[1]-my_values[1]) > margem) ) and  iteracao < iteracoes_maximas ):
        old_values2 = my_values
        hessiana = gradiente_2(funcao, my_values[0], my_values[1])
        d = multiplica_matrix_22_21(hessiana,gradiente(my_values[0],my_values[1],funcao))
        d = (-d[0][0], -d[1][0])
        t = secao_aurea(funcao, epslon, ro, d)
        old_values = my_values
        my_values = (my_values[0] + d[0]*t, my_values[1] + d[1]*t)

        # Inincio BFGS
        p  = [my_values[0]-old_values[0],my_values[1]-old_values[1]]
        grad_delt_x = gradiente(my_values[0],my_values[1],funcao)[0] - gradiente(old_values[0],old_values[1],funcao)[0]
        grad_delt_y = gradiente(my_values[0],my_values[1],funcao)[1] - gradiente(old_values[0],old_values[1],funcao)[1]
        q = (grad_delt_x,grad_delt_y)

        a = hessiana[0][0]
        b = hessiana[0][1]
        c = hessiana[1][0]
        d = hessiana[1][1]

        H = [[(2*p[0]*q[0]*a) + (p[0]*q[1]*(c+b)), (p[0]*q[0]*b) + (p[1]*q[1]*d) + (p[0]*q[1]*d) + (p[1]*q[0]*a)],
            [(p[1]*q[1]*c) + (p[0]*q[0]*c) + (p[0]*q[1]*d) + (p[1]*q[0]*a), 2*p[1]*q[1]*d + p[1]*q[0]*(b+c)]]

        iteracao += 1
        return (my_values,iteracao,old_values[0]-my_values[0],my_function(funcao,my_values[0],my_values[1]))

#my_values = (0.5,0.5)
#d = (-0.004000000000000001, -2.9160000000000004)
#print(metodo_gradiente(3,0.1,0.9))
#print(d[0])

print(metodo_gradiente(1,0.1,0.1))