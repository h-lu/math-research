#!/usr/bin/env python3
"""Arb validated integration of the actual matched harvesting first return.

The branch ODE uses dt/ds=x>0. Every step is enclosed by a Picard tube and
Taylor's theorem: degree N coefficients at the incoming interval plus the
degree N+1 coefficient evaluated on the whole invariant tube. Event-free
steps require interval exclusion or strict monotonicity; crossing events
require opposite endpoint signs and a strict derivative. Interval Newton
isolates their times. No floating tail estimate is used as an error bound.
"""
import argparse
from decimal import Decimal, InvalidOperation
import json
import re
import time
from pathlib import Path
from flint import arb, ctx

NAMES=('beta','d','h','u','v','L','U')

def decimal_literal(value,name):
    """Require one finite exact decimal literal, never a float or Arb ball."""
    if not isinstance(value,str) or re.fullmatch(r'[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?',value) is None:
        raise ValueError(name+' must be a finite exact decimal string')
    try:
        number=Decimal(value)
    except InvalidOperation as error:
        raise ValueError(name+' must be a finite exact decimal string') from error
    if not number.is_finite():
        raise ValueError(name+' must be a finite exact decimal string')
    return number

def validate_settings(dps,order,max_step,max_steps=50000):
    """Reject unsupported settings before changing precision or integrating.

    At least 30 digits keeps the event tolerance 10**(-dps+15) small;
    order two is the supported minimum for the quadratic-extremum start.
    These are API support limits, not a replacement for enclosure checks.
    """
    for name,value,minimum in [('dps',dps,30),('order',order,2),
                               ('max_steps',max_steps,1)]:
        if isinstance(value,bool) or not isinstance(value,int) or value<minimum:
            raise ValueError(name+' must be an integer >= '+str(minimum))
    if decimal_literal(max_step,'max_step')<=0:
        raise ValueError('max_step must be a strictly positive finite exact decimal string')

def sign(x):
    return 1 if x>0 else -1 if x<0 else 0

def conv(a,b,n):
    return sum((a[k]*b[n-k] for k in range(n+1)),arb(0))

def settings(p,region):
    if region=='left': return arb(0),arb(0)
    if region=='middle': return p['u'],p['u']*p['L']
    return p['u']+p['v'],p['u']*p['L']+p['v']*p['U']

def rhs(z,p,region):
    x,y=z[:2]; e,c=settings(p,region)
    return [x*((1-e)*x-x*x-x*y+c),
            p['d']*x*y-p['beta']*y*y-p['h']*x,
            x,(1+p['d']-e)*x-2*x*x-x*y-2*p['beta']*y]

def field_F(z,p,region):
    x,y=z[:2]; e,c=settings(p,region)
    return x*(1-e-x-y)+c

def field_G(z,p):
    x,y=z[:2]
    return p['d']*y-p['beta']*y*y/x-p['h']

def coeffs(z,p,region,nmax):
    e,c=settings(p,region)
    x,y,t,j=([a] for a in z)
    xx=[]; xy=[]; yy=[]
    for n in range(nmax):
        xx.append(conv(x,x,n)); xy.append(conv(x,y,n)); yy.append(conv(y,y,n))
        x.append(((1-e)*xx[n]-conv(xx,x,n)-conv(xx,y,n)+c*x[n])/(n+1))
        y.append((p['d']*xy[n]-p['beta']*yy[n]-p['h']*x[n])/(n+1))
        t.append(x[n]/(n+1))
        j.append(((1+p['d']-e)*x[n]-2*xx[n]-xy[n]-2*p['beta']*y[n])/(n+1))
    return [x,y,t,j]

def polyval(c,s):
    out=c[-1]
    for a in reversed(c[:-1]): out=out*s+a
    return out

class Step:
    def __init__(self,z,p,region,h,order):
        self.z=z; self.p=p; self.region=region; self.h=h; self.order=order
        # Choose an exact dyadic H first, then test a time ball containing
        # [0,H]. This remains valid when h is a decimal ball and upper()
        # itself rounds the old endpoint slightly upward.
        self.valid_time_upper=arb(h/2,h/2).upper()
        if not self.valid_time_upper.is_exact() or not self.valid_time_upper>0:
            raise ArithmeticError('certified time bound is not an exact positive point')
        interval_h=arb(self.valid_time_upper/2,self.valid_time_upper/2)
        f0=rhs(z,p,region)
        # At an extremum x'=0 but x moves quadratically. A quadratic seed
        # prevents the invariant-tube iteration from starting at zero width.
        radii=[2*h*abs(a).upper()+arb('.01')*h*h+arb(10)**(-ctx.dps+5) for a in f0]
        for _ in range(30):
            tube=[arb(a.mid(),a.rad()+r) for a,r in zip(z,radii)]
            images=[a+interval_h*f for a,f in zip(z,rhs(tube,p,region))]
            if all(a.contains_interior(b) for a,b in zip(tube,images)):
                break
            radii=[r*2 for r in radii]
        else: raise ArithmeticError('Picard tube not established')
        if not tube[0]>0 or not tube[1]>0:
            raise ArithmeticError('Picard tube leaves positive quadrant')
        self.tube=tube
        self.c=coeffs([a.mid() for a in z],p,region,order)
        self.remainder=[a[-1] for a in coeffs(tube,p,region,order+1)]
        x,y=tube[:2]; e,c=settings(p,region)
        # Componentwise mean-value comparison removes repeated algebraic
        # dependency on the incoming interval in the Taylor coefficients.
        self.M=[[abs(2*(1-e)*x-3*x*x-2*x*y+c).upper(),abs(x*x).upper()],
                [abs(p['d']*y-p['h']).upper(),abs(p['d']*x-2*p['beta']*y).upper()]]
        self.norm=max(sum(row,arb(0)).upper() for row in self.M)
        self.Mz=[abs(1+p['d']-e-4*x-y).upper(),abs(x+2*p['beta']).upper()]

    def incoming_error(self,s):
        # s may be an event-time ball: all times belong to [0,s.upper()].
        H=abs(s).upper(); A=[[H*a for a in row] for row in self.M]
        r=[a.rad() for a in self.z[:2]]; total=r[:]; term=r[:]
        for k in range(1,25):
            term=[sum((A[i][j]*term[j] for j in range(2)),arb(0))/k for i in range(2)]
            total=[a+b for a,b in zip(total,term)]
        nH=self.norm*H
        tail=nH.exp()*nH**25/arb.fac_ui(25)*max(r)
        total=[(a+tail).upper() for a in total]
        growth=max(r)*nH.exp()
        return [*total,(self.z[2].rad()+H*growth).upper(),
                (self.z[3].rad()+H*sum(self.Mz,arb(0))*growth).upper()]

    def at(self,s):
        # Contract: enclose the actual flow for s in S intersect [0,H], where
        # H=valid_time_upper. Extra times in an overestimated S are harmless
        # polynomial over-evaluation; the remainder is claimed only on [0,H].
        err=self.incoming_error(s)
        return [polyval(c,s)+r*s**(self.order+1)+arb(0,e)
                for c,r,e in zip(self.c,self.remainder,err)]

def event(z,p,region,name):
    if name=='extremum': return field_F(z,p,region)
    return z[0]-p['L' if name.startswith('L') else 'U']

def event_derivative(z,p,region,name):
    v=rhs(z,p,region)
    if name!='extremum': return v[0]
    x,y=z[:2]; e,_=settings(p,region)
    return (1-e-2*x-y)*v[0]-x*v[1]

def event_names(region):
    return (['L-up'] if region=='left' else ['U-down'] if region=='right'
            else ['L-down','U-up'])+['extremum']

def enclose_event(step,name,tol):
    h=step.h
    I=arb(h/2,h/2)
    if not I.contains(0): raise ArithmeticError('initial event-time ball omits zero')
    for _ in range(30):
        mid=I.mid()
        # Arb.intersection encloses an intersection and can slightly overhang
        # an earlier ball. The Newton evaluation center must stay in the
        # certified Picard time domain for the mean-value argument to apply.
        if not 0<=mid<=step.valid_time_upper:
            raise ArithmeticError('event Newton midpoint outside certified time domain')
        zmid=step.at(mid)
        derivative=event_derivative(step.at(I),step.p,step.region,name)
        if not sign(derivative):
            derivative=event_derivative(step.tube,step.p,step.region,name)
        if not sign(derivative): raise ArithmeticError('event derivative unresolved')
        newton=mid-event(zmid,step.p,step.region,name)/derivative
        J=I.intersection(newton)
        if J.rad()>=I.rad()*arb('.99') or J.rad()<tol:
            return J,step.at(J)
        I=J
    raise ArithmeticError('event interval Newton did not settle')

def ball_string(x,dps):
    return x.str(dps,more=True)

def validate_return(parameters,q,section='minimum',dps=70,order=20,max_step='.5',
                    q_radius='0',max_steps=50000):
    """Enclose a return for an exact decimal point or an input ball.

    A successful result certifies existence/uniqueness of the listed
    chronological return for every input point in the ball. It reports
    interval displacement and the exact section derivative formula.
    """
    validate_settings(dps,order,max_step,max_steps)
    if not isinstance(section,str) or section not in ('minimum','maximum'):
        raise ValueError('section must be minimum or maximum')
    if not isinstance(parameters,dict) or set(parameters)!=set(NAMES):
        raise ValueError('parameters must be a dictionary of seven exact decimal strings')
    for name in NAMES: decimal_literal(parameters[name],name)
    decimal_literal(q,'q')
    if decimal_literal(q_radius,'q_radius')<0: raise ValueError('input radius must be nonnegative')
    clock=time.monotonic(); ctx.dps=dps
    p={k:arb(parameters[k]) for k in NAMES}
    if not all(p[k]>0 for k in NAMES) or not 0<p['L']<p['U']<1:
        raise ValueError('physical parameter domain')
    x0=arb(arb(q),arb(q_radius))
    if x0<p['L']: region='left'
    elif x0>p['U']: region='right'
    elif p['L']<=x0 and x0<=p['U']: region='middle'
    else: raise ValueError('input section ball intersects a seam')
    e,c=settings(p,region)
    y0=1-e-x0+c/x0
    initial=[x0,y0,arb(0),arb(0)]
    G0=field_G(initial,p)
    desired_G=-1 if section=='minimum' else 1
    if sign(G0)!=desired_G: raise ValueError('input is not the requested nondegenerate extremum')
    state=initial; departed={'extremum':-desired_G}; extrema=0; logs=[]
    min_positive=arb(1); max_width=arb(0); rejected=0
    base_h=arb(max_step); tol=arb(10)**(-dps+15)
    for count in range(max_steps):
        h=base_h
        while True:
            try:
                st=Step(state,p,region,h,order)
                ending=st.at(h)
                crossings=[]
                for name in event_names(region):
                    ht=event(st.tube,p,region,name)
                    if sign(ht): continue
                    der=event_derivative(st.tube,p,region,name)
                    a=event(state,p,region,name); b=event(ending,p,region,name)
                    if name in departed:
                        if sign(der)!=departed[name]: raise ArithmeticError('departure not monotone in established direction')
                        # Exact event constraint imposed at the previous reset.
                        # Strict monotonicity excludes a second zero in (0,h].
                        continue
                    if sign(a) and sign(b) and sign(a)==sign(b) and sign(der): continue
                    if sign(a) and sign(b) and sign(a)!=sign(b) and sign(der):
                        crossings.append(name); continue
                    raise ArithmeticError('event exclusion unresolved: '+name+
                                          '; count='+str(count)+'; state='+str(state)+
                                          '; a='+str(a)+'; b='+str(b)+'; der='+str(der))
                if len(crossings)>1: raise ArithmeticError('multiple events in one step')
                break
            except ArithmeticError:
                rejected+=1; h/=2
                if h<arb('1e-14'): raise
        min_positive=min_positive.min(st.tube[0].lower()).min(st.tube[1].lower())
        max_width=max_width.max(max(a.rad() for a in ending))
        if not crossings:
            state=ending
            departed={name:direction for name,direction in departed.items()
                      if not sign(event(state,p,region,name))}
            continue
        name=crossings[0]
        root,z=enclose_event(st,name,tol)
        if not root>0: raise ArithmeticError('event root not strictly after step start')
        der=event_derivative(st.tube,p,region,name)
        if name=='extremum':
            extrema+=1
            G=field_G(z,p)
            expected=-desired_G if extrema==1 else desired_G
            if sign(G)!=expected: raise ArithmeticError('wrong extremum order')
            # The exact event is F=0. Project its interval enclosure onto this
            # algebraic nullcline without discarding any exact event point.
            e,c=settings(p,region); z[1]=1-e-z[0]+c/z[0]
            kind='maximum' if expected>0 else 'minimum'
            logs.append({'event':kind,'region':region,
                         'time':ball_string(z[2],dps),'x':ball_string(z[0],dps),
                         'y':ball_string(z[1],dps),'transversality':ball_string(der,dps)})
            if extrema==2:
                derivative=z[3].exp()*G0/field_G(z,p)
                D=z[0]-x0
                return {'status':'VALIDATED actual first return with Arb interval Taylor integration',
                        'parameters':parameters,'q':q,'q_radius':q_radius,'section':section,
                        'dps':dps,'order':order,'max_step':max_step,
                        'P':ball_string(z[0],dps),'displacement':ball_string(D,dps),
                        'displacement_sign':sign(D),
                        'P_prime':ball_string(derivative,dps),
                        'D_prime':ball_string(derivative-1,dps),
                        'D_prime_sign':sign(derivative-1),
                        'physical_period':ball_string(z[2],dps),
                        'positive_tube_margin':ball_string(min_positive,dps),
                        'largest_pre_event_state_radius':ball_string(max_width,dps),
                        'accepted_steps':count+1,'step_halvings':rejected,
                        'events':logs,'elapsed_seconds':time.monotonic()-clock,
                        'claims':['Picard tube and Taylor remainder on every smooth segment.',
                                  'Strict interval exclusion/monotonicity of all enabled events.',
                                  'Unique transverse event roots and actual chronological first return.',
                                  'Whole-tube positivity; all physical matching is exact by model definition.'],
                        'scope':'Point/ball return enclosure; root brackets or global root count require separate checks.'}
            state=z; departed={'extremum':-expected}
        else:
            if sign(der)!=(1 if name.endswith('up') else -1):
                raise ArithmeticError('wrong seam direction')
            z[0]=p['L' if name.startswith('L') else 'U']
            before=region
            region={'L-up':'middle','L-down':'left','U-up':'right','U-down':'middle'}[name]
            logs.append({'event':name,'from':before,'to':region,
                         'time':ball_string(z[2],dps),'x':ball_string(z[0],dps),
                         'y':ball_string(z[1],dps),'transversality':ball_string(der,dps)})
            state=z
            # The newly enabled seam has the same zero with opposite name.
            counterpart={'L-up':'L-down','L-down':'L-up','U-up':'U-down','U-down':'U-up'}[name]
            departed={counterpart:sign(der)}
        if state[2]>300: raise RuntimeError('physical time bound exceeded')
    raise RuntimeError('validated step limit exceeded')

def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('--params',required=True)
    a.add_argument('--q',required=True)
    a.add_argument('--q-radius',default='0')
    a.add_argument('--section',choices=['minimum','maximum'],default='minimum')
    a.add_argument('--dps',type=int,default=70,help='Integer >= 30')
    a.add_argument('--order',type=int,default=20,help='Integer >= 2')
    a.add_argument('--max-step',default='.5',help='Strictly positive finite exact decimal string')
    args=a.parse_args()
    validate_settings(args.dps,args.order,args.max_step)
    source=args.params if args.params.startswith('{') else Path(args.params).read_text()
    print(json.dumps(validate_return(json.loads(source),args.q,args.section,args.dps,args.order,
                                     args.max_step,args.q_radius),indent=2))

if __name__=='__main__': main()
