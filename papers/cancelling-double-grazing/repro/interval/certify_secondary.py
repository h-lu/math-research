#!/usr/bin/env python3
"""Direct fixed-tuple certificate: one negative and five positive q cycles.

This separate witness uses the same harvesting model and the same fixed beta
and u as the main witness. It does not identify the paper's observed eta,d:
the formal secondary target motivated numerical shape fitting only.
"""
import argparse
from decimal import Decimal,localcontext
import json
from pathlib import Path
import platform
import time
import flint
from validated_return import validate_return, validate_settings

CENTERS=[
 '-0.00001706400935651388556564965020739642536671039480734967462752056751636246',
 '0.00000000000176897871053258501269121966725794107122948330929767569622599427632113',
 '0.00000000001320383685043057023513164103967552807810928575633677745603942016869857',
 '0.00000002249999999999999999999999999998381540782246135467363926438981465702725',
 '0.0000002025000000000000000000000000000126889345860275331106281694169359403343',
 '0.000210982484000851405019370526095347199859928350259824810116963158712512',
]
RADIUS='0.000000000000000000000000000000000001'

def difference(a,b):
    with localcontext() as context:
        context.prec=180
        return format(Decimal(a)-Decimal(b),'f')

def sum_decimal(a,b):
    with localcontext() as context:
        context.prec=180
        return format(Decimal(a)+Decimal(b),'f')

def certify(parameters,dps=120,order=60,max_step='.5',progress=False):
    validate_settings(dps,order,max_step)
    started=time.monotonic()
    brackets=[[difference(c,RADIUS),sum_decimal(c,RADIUS)] for c in CENTERS]
    if not(Decimal(brackets[0][1])<0 and all(Decimal(a)>0 for a,b in brackets[1:])):
        raise ArithmeticError('one-negative five-positive classification failed')
    if not all(Decimal(brackets[i][1])<Decimal(brackets[i+1][0]) for i in range(5)):
        raise ArithmeticError('isolators overlap')
    results=[]
    for index,(center,bracket) in enumerate(zip(CENTERS,brackets),1):
        endpoints=[]
        for q in bracket:
            r=validate_return(parameters,difference(parameters['L'],q),section='minimum',
                              dps=dps,order=order,max_step=max_step)
            r['outward_q']=q
            r['outward_displacement_sign']=-r['displacement_sign']
            endpoints.append(r)
        if endpoints[0]['displacement_sign']*endpoints[1]['displacement_sign']!=-1:
            raise ArithmeticError('opposite endpoint signs not proved: '+str(index))
        ball=validate_return(parameters,difference(parameters['L'],center),section='minimum',
                             dps=dps,order=order,max_step=max_step,q_radius=RADIUS)
        if not ball['D_prime_sign']:
            raise ArithmeticError('whole-isolator derivative exclusion failed: '+str(index))
        expected=(['maximum','minimum'] if index==1 else
                  ['L-up','maximum','L-down','minimum'] if index in (2,3) else
                  ['L-up','U-up','maximum','U-down','L-down','minimum'])
        if any([e['event'] for e in r['events']]!=expected for r in [*endpoints,ball]):
            raise ArithmeticError('unexpected chronology: '+str(index))
        results.append({'index':index,'outward_q_bracket':bracket,'radius':RADIUS,
                        'endpoint_enclosures':endpoints,'whole_bracket_return':ball,
                        'stability':'repelling' if ball['D_prime_sign']>0 else 'attracting',
                        'conclusion':'Exactly one simple fixed point in this isolator.'})
        if progress: print('Certified secondary isolator',index,'Dprime',ball['D_prime'],flush=True)
    return {'status':'VALIDATED_AT_LEAST_SIX_SIMPLE_ACTUAL_PERIODIC_ORBITS_ONE_NEGATIVE_FIVE_POSITIVE',
            'parameters':parameters,'parameter_interpretation':'Exact decimal rational values.',
            'coordinate':'q=L-x_minimum; Dq=-(Px-x). Signs refer to this explicitly fixed physical coordinate.',
            'integration':{'arithmetic':'python-flint Arb','dps':dps,'order':order,
                           'max_rescaled_time_step':max_step,'python_version':platform.python_version(),
                           'python_flint_version':flint.__version__,'flint_version':flint.__FLINT_VERSION__},
            'isolators':results,'certified_count_lower_bound':6,
            'certified_negative_q_isolators':1,'certified_positive_q_isolators':5,
            'exact_total_count_certified':False,
            'proof':['Each disjoint rational isolator has opposite rigorous endpoint signs.',
                     'Each whole isolator has one common chronological return with positive tubes and transverse events.',
                     'Pprime-1 excludes zero on the whole isolator, proving one simple fixed point.',
                     'One isolator lies at q<0 and five at q>0, giving six distinct actual periodic orbits.'],
            'limitations':['No exact total upper bound or whole-interval exclusion.',
                           'The formal eta=2/(27B^2), d=1/(8B) target motivated shape fitting; no observed-coordinate identity is proved.',
                           'No exact neutral or cancelling reference is asserted.'],
            'elapsed_seconds':time.monotonic()-started}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--params',default=str(Path(__file__).with_name('secondary-candidate.json')))
    parser.add_argument('--output',default=None)
    parser.add_argument('--dps',type=int,default=120,help='Integer >= 30')
    parser.add_argument('--order',type=int,default=60,help='Integer >= 2')
    parser.add_argument('--max-step',default='.5',help='Strictly positive finite exact decimal string')
    parser.add_argument('--progress',action='store_true')
    args=parser.parse_args()
    validate_settings(args.dps,args.order,args.max_step)
    parameters=json.loads(Path(args.params).read_text())
    report=certify(parameters,args.dps,args.order,args.max_step,args.progress)
    payload=json.dumps(report,indent=2)+'\n'
    if args.output: Path(args.output).write_text(payload,encoding='utf-8')
    else: print(payload,end='')

if __name__=='__main__':main()
