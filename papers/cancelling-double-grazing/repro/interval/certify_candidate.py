#!/usr/bin/env python3
"""Fail-closed certificate of at least six simple actual harvesting cycles.

All seven physical parameters and bracket endpoints are exact decimal
rationals. This program certifies the actual piecewise model, without a
neutral-reference or observed-normal-form assertion. It does not prove an
upper bound on the number of cycles outside these six disjoint isolators.
"""
import argparse
from decimal import Decimal, localcontext
import json
from pathlib import Path
import platform
import time
import flint
from validated_return import validate_return, validate_settings

PARAMETERS={
    'beta':'0.028676736587811624',
    'd':'0.4374783057744852',
    'h':'0.25790748326214347',
    'u':'0.01',
    'v':'0.003857453910355628',
    'L':'0.170270029350295',
    'U':'0.2549289376256445',
}
CENTERS=[
    '-0.004256180956417886653696195992805154317281',
    '-0.00001517561150190900562888108319718409583430',
    '0.00002500000305083428568701180976919153299133',
    '0.00009999999919525798963346717645775552106603',
    '0.0002249999661036241957609848382948107187517',
    '0.0004000000224609623026298208486068033867455',
]
RADIUS='0.000000000000000000000001'

def decimal_arithmetic(a,b,subtract=False):
    with localcontext() as context:
        context.prec=120
        return format(Decimal(a)-Decimal(b) if subtract else Decimal(a)+Decimal(b),'f')

def certify(dps=100,order=48,max_step='.5',progress=False):
    validate_settings(dps,order,max_step)
    started=time.monotonic()
    brackets=[[decimal_arithmetic(c,RADIUS,True),decimal_arithmetic(c,RADIUS)]
              for c in CENTERS]
    if not all(Decimal(a)<Decimal(b) for a,b in brackets):
        raise ArithmeticError('a bracket has nonpositive width')
    if not all(Decimal(brackets[i][1])<Decimal(brackets[i+1][0]) for i in range(5)):
        raise ArithmeticError('root isolators overlap')
    results=[]
    for index,(center,bracket) in enumerate(zip(CENTERS,brackets),1):
        endpoint_results=[]
        for q in bracket:
            x=decimal_arithmetic(PARAMETERS['L'],q,True)
            r=validate_return(PARAMETERS,x,section='minimum',dps=dps,
                              order=order,max_step=max_step)
            r['outward_q']=q
            r['outward_displacement_sign']=-r['displacement_sign']
            endpoint_results.append(r)
        if endpoint_results[0]['displacement_sign']*endpoint_results[1]['displacement_sign']!=-1:
            raise ArithmeticError('opposite rigorous endpoint signs not proved for bracket '+str(index))
        xcenter=decimal_arithmetic(PARAMETERS['L'],center,True)
        derivative_ball=validate_return(PARAMETERS,xcenter,section='minimum',dps=dps,
                                        order=order,max_step=max_step,q_radius=RADIUS)
        if not derivative_ball['D_prime_sign']:
            raise ArithmeticError('whole-isolator derivative exclusion not proved for bracket '+str(index))
        expected_events=(['maximum','minimum'] if index<=2 else
                         ['L-up','U-up','maximum','U-down','L-down','minimum'])
        if any([e['event'] for e in r['events']]!=expected_events
               for r in [*endpoint_results,derivative_ball]):
            raise ArithmeticError('chronology differs on bracket '+str(index))
        item={'index':index,'outward_q_bracket':bracket,'radius':RADIUS,
              'endpoint_enclosures':endpoint_results,'whole_bracket_return':derivative_ball,
              'conclusion':'Exactly one simple fixed point in this isolator.',
              'stability':'repelling' if derivative_ball['D_prime_sign']>0 else 'attracting'}
        results.append(item)
        if progress:
            print('Certified isolator',index,'Dprime',derivative_ball['D_prime'],flush=True)
    return {'status':'VALIDATED_AT_LEAST_SIX_SIMPLE_ACTUAL_PERIODIC_ORBITS',
            'parameters':PARAMETERS,'parameter_interpretation':'Exact decimal rational values.',
            'coordinate':'q=L-x at the actual physical prey-minimum nullcline; Dq=-(Px-x).',
            'integration':{'arithmetic':'python-flint Arb','dps':dps,'order':order,
                           'max_rescaled_time_step':max_step,
                           'python_version':platform.python_version(),
                           'python_flint_version':flint.__version__,
                           'flint_version':flint.__FLINT_VERSION__},
            'isolators':results,
            'certified_count_lower_bound':6,'exact_total_count_certified':False,
            'proof':['Each endpoint return has strict opposite displacement sign in its isolator.',
                     'The whole isolator has a common positive-domain chronological first return with transverse events.',
                     'Its return derivative minus one excludes zero throughout the isolator.',
                     'Continuity and the intermediate value theorem yield one fixed point; derivative exclusion gives uniqueness and simplicity.',
                     'The six disjoint minimum-section isolators give six distinct periodic orbits.'],
            'limitations':['No whole-interval exclusion or verified analytic upper bound: the total count is at least six.',
                           'No claim that the rounded smooth-fold reference is exactly neutral or cancelling.',
                           'No observed normal-form coordinates are identified.'],
            'elapsed_seconds':time.monotonic()-started}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',default=None)
    parser.add_argument('--dps',type=int,default=100,help='Integer >= 30')
    parser.add_argument('--order',type=int,default=48,help='Integer >= 2')
    parser.add_argument('--max-step',default='.5',help='Strictly positive finite exact decimal string')
    parser.add_argument('--progress',action='store_true')
    args=parser.parse_args()
    # Failure raises and never emits the certified status.
    report=certify(args.dps,args.order,args.max_step,args.progress)
    payload=json.dumps(report,indent=2)+'\n'
    if args.output: Path(args.output).write_text(payload,encoding='utf-8')
    else: print(payload,end='')

if __name__=='__main__':main()
