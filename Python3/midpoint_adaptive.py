#! /usr/bin/env python3
#
def lotka_deriv ( t, y ):

#*****************************************************************************80
#
## lotka_deriv() evaluates the right hand side of a Lotka-Volterra ODE.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    29 May 2022
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real T: the time
#
#    real Y(1,2): the current solution values.
#
#  Output:
#
#    real DYDT(2,1): the derivative values.
#
  import numpy as np

  u = y[0]
  v = y[1]

  dudt =  - 2.0 * u + u * v
  dvdt =          v - u * v

  dydt = np.array ( [ dudt, dvdt ] )

  return dydt

def lotka_test ( ):

#*****************************************************************************80
#
## lotka_test(): midpoint_adaptive() solves lotka_ode().
#
#  Discussion:
#
#    The physical system under consideration is a pair of animal populations.
#
#    The PREY reproduce rapidly for each animal alive at the beginning of the
#    year, two more will be born by the end of the year.  The prey do not have
#    a natural death rate instead, they only die by being eaten by the predator.
#    Every prey animal has 1 chance in 1000 of being eaten in a given year by
#    a given predator.
#
#    The PREDATORS only die of starvation, but this happens very quickly.
#    If unfed, a predator will tend to starve in about 1/10 of a year.
#    On the other hand, the predator reproduction rate is dependent on
#    eating prey, and the chances of this depend on the number of available prey.
#
#    The resulting differential equations can be written:
#
#      PREY(0) = 5000         
#      PRED(0) =  100
#
#      d PREY / dT =    2 * PREY(T) - 0.001 * PREY(T) * PRED(T)
#      d PRED / dT = - 10 * PRED(T) + 0.002 * PREY(T) * PRED(T)
#
#    Here, the initial values (5000,100) are a somewhat arbitrary starting point.
#
#    The pair of ordinary differential equations that result have an interesting
#    behavior.  For certain choices of the interaction coefficients (such as
#    those given here), the populations of predator and prey will tend to 
#    a periodic oscillation.  The two populations will be out of phase the number
#    of prey will rise, then after a delay, the predators will rise as the prey
#    begins to fall, causing the predator population to crash again.
#
#    There is a conserved quantity, which here would be:
#      E(r,f) = 0.002 r + 0.001 f - 10 ln(r) - 2 ln(f)
# 
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    29 May 2022
#
#  Author:
#
#    John Burkardt
#
#  Reference:
#
#    George Lindfield, John Penny,
#    Numerical Methods Using MATLAB,
#    Second Edition,
#    Prentice Hall, 1999,
#    ISBN: 0-13-012641-1,
#    LC: QA297.P45.
#
#  Input:
#
#    real TSPAN = [ T0, TMAX ], the initial and final times.
#    A reasonable value might be [ 0, 5 ].
#
#    real Y0 = [ PREY, PRED ], the initial number of prey and predators.
#    A reasonable value might be [ 5000, 100 ].
#
  import matplotlib.pyplot as plt
  import numpy as np

  print ( '' )
  print ( 'lotka_test():' )
  print ( '  midpoint_adaptive() solves lotka_ode().' )
  print ( '  A pair of ordinary differential equations for a population' )
  print ( '  of predators and prey are solved using midpoint().' )
  print ( '' )
  print ( '  The exact solution shows periodic behavior, with a fixed' )
  print ( '  period and amplitude.' )

  f = lotka_deriv
  a = 0.0
  b = 100.0
  ya = np.array ( [ 4.0, 2.0 ] )
  tau0 = 2.5e-2
  reltol = 1.0e-3
  abstol = 1.0e-3

  t, y = mad_compute ( f, a, b, ya, tau0, reltol, abstol )

  mad_stats ( t, 'lotka' )
#
#  Time plot.
#
  plt.plot ( t, y[:,0], 'r-', linewidth = 2 )
  plt.plot ( t, y[:,1], 'g-', linewidth = 2 )
  plt.grid ( True )
  plt.xlabel ( '<--- Prey --->' )
  plt.ylabel ( '<--- Predators --->' )
  plt.title ( 'midpoint(): lotka time plot' )
  plt.legend ( [ 'prey', 'predator' ] )
  filename = 'lotka_time.png'
  plt.savefig ( filename )
  print ( '  Graphics saved as "%s"' % ( filename ) )
  plt.show ( block = False )
  plt.close ( )
#
#  Phase plot.
#
  plt.plot ( y[:,0], y[:,1], 'r-', linewidth = 2 )
  plt.grid ( True )
  plt.xlabel ( '<--- Prey --->' )
  plt.ylabel ( '<--- Predators --->' )
  plt.title ( 'midpoint(): lotka phase plot' )
  filename = 'lotka_phase.png'
  plt.savefig ( filename )
  print ( '  Graphics saved as "%s"' % ( filename ) )
  plt.show ( block = False )
  plt.close ( )

  return

def mad_compute ( f, a, b, ya, tau0, reltol, abstol ):

#*****************************************************************************80
#
## mad_compute() estimates the solution of an ordinary differential equation.
#
#  Discussion:
#
#    mad_compute() uses an adaptive time stepping algorithm for the
#    implicit midpoint method, applied to a system of ordinary differential
#    equations.
#
#    The implicit midpoint method has local truncation error O(h^3).
# 
#    The time step is adaptively controlled with respect to relative and
#    absolute tolerances applied to the estimated local truncation error (LTE).
#    
#    The adaptivity with respect to relative and absolute error tolerances
#    is described on page 168 in the Hairer reference.
#
#    The function can be invoked by a command like:
#
#      [ t, y ] = mad_compute ( @ lotka_deriv, 0.0, 250.0, [6.0,2.0], ...
#         25.e-3, 1.e-4, 1.e-4 )
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    29 May 2022
#
#  Author:
#
#    John Burkardt, Catalin Trenchea
#
#  Reference:
#
#    Burkardt, John Trenchea, Catalin,
#    Refactorization of the midpoint rule,
#    Appl. Math. Lett. 107 (2020), 106438.
#
#    Hairer, E. Nørsett, S. P. Wanner, G. 
#    Solving ordinary differential equations, I. Nonstiff problems, 
#    Springer Series in Computational Mathematics, Number 8. 
#    Springer-Verlag, Berlin, 1987.
#
#  Input:
#
#    function handle f: the function that evaluates the right hand side 
#    of the differential equation, of the form 
#      dydt = f ( t, y )
#    The user can define the function in an M file, and input f as the
#    name of this file, preceded by an @ sign.  An alternative for
#    experienced users, and simple right hand sides, is to use an
#    anonymous function, in which case f might be passed as
#    '@(x)2.0*sin(x)', for example.
#
#    real a, b: the interval of integration.
#
#    real ya(m): the (row) vector of initial conditions at the starting time. 
#
#    real tau0: the initial timestep to take.
#
#    real reltol, abstol: the tolerances for the local truncation error.
#
#  Output:
#
#    real t(n+1), y(n+1,m): the sequence of times and solution estimates.
#
#  Local:
#
#    integer count: the number of stepsizes, whether accepted or rejected.
#    It should be the case that it <= count.
#
#    real kappa: a factor for the stepsize update.
#    0 < kappa <= 1.  The value kappa = 0.85 is suggested.
#
#    real theta: the theta-method parameter.
#    theta = 0.0 for the backward Euler method.
#    theta = 0.5 for the midpoint method.
#    theta = 1.0 for the forward Euler method.
#
  import numpy as np

  kappa = 0.85
  theta = 0.5
  m = len ( ya )

  it = -2
  count = -1
  t = np.zeros ( 0 )
  y = np.zeros ( [ 1, m ] )
  tau = np.zeros ( 0 )

  while ( True ):

    it = it + 1
#
#  it = -1: y0.
#
    if ( it == -1 ):
      t = np.append ( t, a )
      y[0,:] = ya[:]
      tau = np.append ( tau, tau0 )
      count = count + 1

    elif ( b <= t[-1] ):
 
      break
#
#  it = 0, 1: y1, y2
#
    elif ( it <= 1 ):
      t, y = mad_step ( it, t, y, tau, f, theta )
      tau = np.append ( tau, tau0 )
      count = count + 1

    else:
#
#  Several step reductions may be necessary before the LTE test is met.
#
      while ( True ):

        count = count + 1

        t, y = mad_step ( it, t, y, tau, f, theta )
#
#  Estimate local truncation error.
# 
        tnmin, tnmax = mad_lte ( it, y, tau, reltol, abstol )
#
#  If the LTE test was met, we accept T and Y, set the next time step,
#  and can advance to the next time.
#
        if ( tnmin <= 1.0 ):
          factor = kappa * ( 1.0 / tnmax ) ** ( 1.0 / 3.0 )
          factor = min ( factor, 1.5 )
          factor = max ( factor, 0.02 )
          tau = np.append ( tau, factor * tau[it] )
          break

  print ( '' )
  print ( '  mad_compute(): count = ', count )

  return t, y

def mad_lte ( it, y, tau, reltol, abstol ):

#*****************************************************************************80
#
## mad_lte() estimates the local truncation error.
#
#  Discussion:
#
#    This function is called by mad_compute() to estimate the local truncation
#    error incurred during a single step of the implicit midpoint method.
#
#    We assume the step has been taken over an interval which we will 
#    represent as going from (t1,y1) to (t2,y2).  We may regard y2 as the
#    estimated solution of the local ODE:
#      dy/dt = f(t,y), y(t1) = y1.
#    We suppose that a function y*() is the EXACT solution of this system
#    and we define the local truncation error for this step as
#      lte = y*(t2) - y2
#    We wish to accept this step as long as lte is small, that is:
#      lte = y*(t2) - y2 < abstol + reltol * max ( y1, y2 )
#    Of course, to carry out this check, we need to estimate lte or y*(t2).
#    This is done using the Milne Device, described in the Milne reference.
#    Essentially, lte is estimated using a weighted difference of solution
#    estimates y2 and an estimate formed by an Adams-Bashforth AB2-like 
#    method.  See page 5 of the Burkardt reference for details.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    29 May 2022
#
#  Author:
#
#    John Burkardt, Catalin Trenchea
#
#  Reference:
#
#    Burkardt, John Trenchea, Catalin,
#    Refactorization of the midpoint rule,
#    Appl. Math. Lett. 107 (2020), 106438.
#
#    Milne, W. E.
#    Numerical Integration of Ordinary Differential Equations,
#    Amer. Math. Monthly,
#    33 (1926), no. 9, 455–460.
#
#  Input:
#
#    integer it: the index of the most recently accepted step.
#
#    real y(it+1,m): the solution vector.
#
#    real tau[it]: the stepsize vector.
#
#    real reltol, abstol: the tolerances for the local truncation error.
#
#  Output:
#
#    real tnmin, tnmax: the minimum and maximum entries of an lte
#    error estimator.
#
  import numpy as np
#
#  Evaluate the AB2-like solution.
#
  c1 = ( tau[it] + tau[it-1] ) * ( tau[it] + tau[it-1] + tau[it-2] ) \
    / ( tau[it-1] * ( tau[it-1] + tau[it-2] ) )

  c2 = - tau[it] * ( tau[it] + tau[it-1] + tau[it-2] ) \
    / ( tau[it-1] * tau[it-2] )

  c3 = tau[it] * ( tau[it] + tau[it-1] ) \
     / ( tau[it-2] * ( tau[it-1] + tau[it-2] ) )

  uab2 = c1 * y[it,:] + c2 * y[it-1,:] + c3 * y[it-2,:]
#
#  Evaluate R, the variable error coefficient in the LTE.
#
  r = 1.0 / 24.0 + 1.0 / 8.0 * ( 1.0 + tau[it-1] / tau[it] ) \
    * ( 1.0 + 2.0 * tau[it-1] / tau[it] + tau[it-2] / tau[it] )
#
#  Use AB2 estimator to approximate the LTE vector.
#
  tn = ( y[it+1,:] - uab2 ) * 24.0 * r / ( 24.0 * r - 1.0 ) \
    / ( abstol + np.abs ( y[it+1,:] ) * reltol )
   
  tnmax = np.max ( np.abs ( tn ) )
  tnmin = np.min ( np.abs ( tn ) )

  return tnmin, tnmax

def mad_stats ( t, label ):

#*****************************************************************************80
#
## mad_stats() reports time and timestep statistics for mad_compute().
#
#  Discussion:
#
#    This function prints some information about the time nodes and time steps,
#    and creates plots that show the variation in timesteps, and in the ratio
#    of successive timesteps, over the integration interval.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    29 May 2022
#
#  Author:
#
#    John Burkardt, Catalin Trenchea
#
#  Input:
#
#    real t(n+1): the sequence of time values.
#
#    string label: a label to prefix the plotfile names.  Default is 'mad'.
#    The files will be 'label_timesteps.png' and 'label_timestep_ratios.png'.
#
  import matplotlib.pyplot as plt
  import numpy as np

  tau = np.diff ( t )
  r = tau[1:] / tau[0:-1]

  print ( '' )
  print ( 'mad_stats():' )
  print ( '' )
  print ( '  Final time interval is [', t[-2], ',', t[-1], ']' )
  print ( '  Number of times t() ', len ( t ) )
  print ( '  Number of timesteps tau() ', len ( tau ) )
  print ( '  Minimum timestep is ', np.min ( tau ) )
  print ( '  Maximum timestep is ', np.max ( tau ) )

  plt.clf ( )
  plt.plot ( t[0:-1], tau, '.-' )
  plt.grid ( True )
  plt.title ( 'timesteps' )
  plt.xlabel ( '<-- t(i) -->' )
  plt.ylabel ( '<-- dt(i) -->' )
  filename = label +'_timesteps.png'
  plt.savefig ( filename )
  print ( '  Graphics saved as "' + filename + '"' )
  plt.show ( block = False )
  plt.close ( )

  plt.clf ( )
  plt.plot ( t[0:-2], r, '.-' )
  plt.grid ( True )
  plt.title ( 'timestep ratios' )
  plt.xlabel ( '<-- t(i) -->' )
  plt.ylabel ( '<-- dt(i+1)/dt(i) -->' )
  filename = label + '_timestep_ratios.png'
  plt.savefig ( filename )
  print ( '  Graphics saved as "' + filename + '"' )
  plt.show ( block = False )
  plt.close ( )

  return

def mad_step ( it, t, y, tau, f, theta ):

#*****************************************************************************80
#
## mad_step extends the (t,y) solution sequence by one more step.
#
#  Discussion:
#
#    Given a sequence of solution values (t,y) for an ODE, this function
#    starts at time t[it], with solution y[it,:] and stepsize tau[it].
#
#    It uses the implicit midpoint method to estimate a solution ym[:] at
#    time tm = t[it] + theta * tau[it], from which it constructs
#    the values t[it+1] and y[it+1,:].
#
#    In other words, it adds one more (t,y) pair to the solution sequence.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    29 May 2022
#
#  Author:
#
#    John Burkardt, Catalin Trenchea
#
#  Reference:
#
#    John Burkardt, Catalin Trenchea,
#    Refactorization of the midpoint rule,
#    Appl. Math. Lett. 107 (2020), 106438.
#
#  Input:
#
#    integer it: the current number of time steps.
#
#    real t[it], y(it,m): the sequence of times and solution estimates.
#
#    real tau[it]: the sequence of time steps.  The final entry is tentative.
#
#    function f: the function that evaluates the right hand side 
#    of the differential equation, of the form 
#      dydt = f ( t, y )
#
#    real theta: the theta-method parameter.
#    theta = 0.0 for the backward Euler method.
#    theta = 0.5 for the midpoint method.
#    theta = 1.0 for the forward Euler method.
#
#  Output:
#
#    real t[it+1], y[it+1,m]: the sequence of times and solution estimates,
#    augmented by one new time step.
#
  from scipy.optimize import fsolve
  import numpy as np

  tm = t[it] + theta * tau[it]
  ym = y[it,:] + theta * tau[it] * f ( tm, y[it,:] )
  ym = fsolve ( midpoint_residual, ym, args = ( f, t[it], y[it,:], tm ) )

  tn = t[it] + tau[it]
  yn = (       1.0 / theta ) * ym \
     + ( 1.0 - 1.0 / theta ) * y[it,:]

  t = np.append ( t, tn )
  y = np.vstack ( ( y, yn ) )

  return t, y

def midpoint_residual ( yh, f, to, yo, th ):

#*****************************************************************************80
#
## midpoint_residual() evaluates the midpoint residual.
#
#  Discussion:
#
#    We are seeking a value YH defined by the implicit equation:
#
#      YH = YO + ( TH - TO ) * F ( TH, YH )
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    29 May 2022
#
#  Author:
#
#    John Burkardt
#
#  Input:
#
#    real yh: the estimated solution value at the midpoint time.
#
#    function f: evaluates the right hand side of the ODE.  
#
#    real to, yo: the old time and solution value.
#
#    real th: the midpoint time.
#
#  Output:
#
#    real value: the midpoint residual.
#
  value = yh - yo - ( th - to ) * f ( th, yh )

  return value

def midpoint_adaptive_test ( ):

#*****************************************************************************80
#
## midpoint_adaptive_test() tests midpoint_adaptive().
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    29 May 2022
#
#  Author:
#
#    John Burkardt
#
  import platform

  print ( '' )
  print ( 'midpoint_adaptive_test():' )
  print ( '  Python version: ' + platform.python_version ( ) )
  print ( '  Test midpoint_adaptive().' )

  lotka_test ( )
#
#  Terminate.
#
  print ( '' )
  print ( 'midpoint_adaptive_test():' )
  print ( '  Normal end of execution.' )

  return

def timestamp ( ):

#*****************************************************************************80
#
## timestamp() prints the date as a timestamp.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license. 
#
#  Modified:
#
#    21 August 2019
#
#  Author:
#
#    John Burkardt
#
  import time

  t = time.time ( )
  print ( time.ctime ( t ) )

  return

if ( __name__ == "__main__" ):
  timestamp ( )
  midpoint_adaptive_test ( )
  timestamp ( )
