"""Reproducible standard-input conversion and native FreeFlux Studio batch run.

The archived stoichiometry is already oriented in the solved direction. Never
invert it again for reactions whose source metadata says 'reverse'.
"""
from __future__ import annotations
import argparse
import copy
import contextlib
import hashlib
import json
import math
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
# Keep the frozen nine-path source inside the repository fixture directory so a
# fresh clone can reproduce the historical conversion without a workstation-
# specific download path.
SOURCE = HERE / 'source' / 'pathways.jsonl'
SOURCE_SHA = 'd7668f5fcbe2c09cdcde3983351fa5fd3dea23bd0af18188c66fdaf309adb0d8'
sys.path.insert(0, str(ROOT))

ALIASES = {
 'bigg:2ddg6p_c':'KDPG', 'bigg:6pgc_c':'SixPG', 'bigg:6pgl_c':'SixPGL',
 'bigg:e4p_c':'E4P', 'bigg:f6p_c':'F6P', 'bigg:fald_c':'Formaldehyde',
 'bigg:g3p_c':'G3P', 'bigg:g6p_c':'G6P', 'bigg:h2o_c':'H2O', 'bigg:h_c':'H',
 'bigg:nadp_c':'NADP', 'bigg:nadph_c':'NADPH', 'bigg:pyr_c':'Pyruvate',
 'bigg:r5p_c':'R5P', 'bigg:ru5p__D_c':'Ru5P', 'bigg:s7p_c':'S7P',
 'bigg:xu5p__D_c':'Xu5P', 'rump:C06019_c':'H6P',
 'bigg:adp_c':'ADP', 'bigg:atp_c':'ATP', 'bigg:dhap_c':'DHAP',
 'bigg:fdp_c':'FBP', 'bigg:lac__D_c':'D_lactate', 'bigg:mthgxl_c':'Methylglyoxal',
 'bigg:nad_c':'NAD', 'bigg:nadh_c':'NADH', 'bigg:pi_c':'Pi',
 'bigg:13dpg_c':'BPG', 'bigg:2pg_c':'PG2', 'bigg:3pg_c':'PG3',
 'bigg:dha_c':'DHA', 'bigg:pep_c':'PEP', 'bigg:fru_c':'Fructose',
 'bigg:glc__D_c':'Glucose', 'bigg:3php_c':'PHP',
 'bigg:glyc__R_c':'Glycerate', 'bigg:hpyr_c':'HPYR', 'bigg:2h3oppan_c':'TSA',
}
DISPLAY = {'SixPG':'6PG', 'SixPGL':'6PGL', 'D_lactate':'D-lactate',
 'BPG':'1,3-BPG', 'PG2':'2PG', 'PG3':'3PG', 'PHP':'3PHP',
 'NAD':'NAD⁺', 'NADP':'NADP⁺', 'H':'H⁺', 'H2O':'H₂O'}
CURRENCY = {'ATP','ADP','NAD','NADH','NADP','NADPH','Pi','H','H2O'}
R_LABEL = {'RUMP_R05338':'Hps', 'RUMP_R05339':'Phi', 'IML_TALA':'TAL',
 'IML_G6PDH2r':'G6PDH', 'IML_LDH_D':'LDH-D', 'IML_GLYOX3':'GLYOX3'}
TITLES = {1:'RuMP-ED',2:'RuMP / methylglyoxal',3:'RuMP / F6PA',
 4:'RuMP / hexose-ED',5:'RuMP / EMP',6:'RuMP / DHAP-DHA',
 7:'RuMP / HPYRRx',8:'RuMP / HPYRRy',9:'RuMP / TSA'}
ED = ['F6P','IML_PGI','G6P','IML_G6PDH2r','SixPGL','IML_PGL','SixPG','IML_EDD','KDPG','IML_EDA','G3P']
HEX = ['F6P','IML_F6PP','Fructose','IML_XYLI2','Glucose','IML_HEX1','G6P','IML_G6PDH2r','SixPGL','IML_PGL','SixPG','IML_EDD','KDPG','IML_EDA','G3P']
EMP = ['F6P','IML_PFK','FBP','IML_FBA','G3P']
F6PA = ['F6P','IML_F6PA','G3P']
REGEN = ['IML_TKT2','E4P','IML_TALA','S7P','IML_TKT1','R5P','IML_RPI']
LOWER = ['G3P','IML_GAPD','BPG','IML_PGK','PG3','IML_PGM','PG2','IML_ENO','PEP','IML_DHAPT','DHAP','IML_TPI']
HP = ['G3P','IML_GAPD','BPG','IML_PGK','PG3','IML_PGCD','PHP','IML_HPYRP','HPYR']
TAIL = ['Glycerate','IML_GLYCK2','PG2','IML_ENO','PEP','IML_DHAPT','DHAP','IML_TPI']

def dump(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')

def point(center, radius, degrees):
    a = math.radians(degrees)
    return [center[0]+radius*math.cos(a), center[1]+radius*math.sin(a)]

def bezier(center, radius, begin, end):
    a, b = math.radians(begin), math.radians(end)
    k = 4/3*math.tan((b-a)/4)
    p, q = point(center, radius, begin), point(center, radius, end)
    return [[p[0]-k*radius*math.sin(a),p[1]+k*radius*math.cos(a)],
            [q[0]+k*radius*math.sin(b),q[1]-k*radius*math.cos(b)]]

def label_offset(name, angle, font=19, radius=9):
    # Account for the full label, including its numeric prefix and pool code.
    width = (len(name)+7)*font*4/3*.53
    height = font*4/3*1.1
    a = math.radians(angle)
    x, y = math.cos(a), math.sin(a)
    distance = min((width/2+radius+12)/max(abs(x),.001),
                   (height/2+radius+12)/max(abs(y),.001))
    return {'label_dx':distance*x, 'label_dy':distance*y-4}

def make_layout(index, pathway, side_pathway_mode='open_chain'):
    reactions = {r['reaction_id']: r for r in pathway['reactions']}
    mids = sorted({ALIASES[m] for r in reactions.values() for m in r['full_stoichiometry']})
    output = ED if index==1 else HEX if index==4 else EMP if index in (2,5,6) else F6PA
    items = ['Ru5P','RUMP_R05338','H6P','RUMP_R05339'] + output + REGEN
    radius = 600 if index==4 else 500
    center = [850,900] if index==4 else [780,820]
    width = 1750 if index==4 else 1600 if index==1 else 3200
    height = 1870 if index==4 else 1710 if index==1 else 1870
    start = 82 - 360*items.index('G3P')/len(items)
    angles = {mid:start+360*i/len(items) for i,mid in enumerate(items)}
    shared_g3p = point(center,radius-175,angles['S7P'])
    shared_xu5p = point(center,150,(angles['IML_TKT2']+angles['IML_TKT1'])/2)
    ru5p = point(center,radius,angles['Ru5P'])
    overrides = {}
    for mid in items[::2]:
        overrides[mid] = label_offset(DISPLAY.get(mid,mid),angles[mid])
    for rid in items[1::2]:
        a=math.radians(angles[rid])
        overrides[rid]={'label_dx':-63*math.cos(a),'label_dy':-63*math.sin(a)-5}
    spec = {
      'enabled':True,'center':center,'radius':radius,'start_angle_deg':start,'items':items,
      'external_metabolites':['Formaldehyde','Pyruvate'] + (['DHAP'] if index in (2,5,6) else []) + (['DHA'] if index in (3,7,8,9) else []),
      'satellite_inset':150,'satellite_spacing':75,
      'shared_metabolites':{'Xu5P':{'position':shared_xu5p,'label_offset':[70,-15]}},
      'shared_instances':{'G3P_regeneration':{'metabolite':'G3P','reactions':['IML_TALA','IML_TKT1'],'position':shared_g3p,'label_offset':[80,-10]}},
      'instance_marker_style':'pool_code','input_arrows':True,
      'side_reactions':[{'reaction':'IML_RPE','source':'Xu5P','target':'Ru5P',
        'control_points':[[center[0]-220,center[1]-80],[ru5p[0]-35,ru5p[1]+120]],'label_offset':[55,-20]}],
      'title':TITLES[index], 'subtitle':f'Candidate {index:02d}',
      'title_position':[center[0]+110,center[1]-135],
      'instance_legend_position':[width/2,height-135], 'view_overrides':overrides,
    }
    layout = {
      'template':'candidate01_demo/empty_template.json','template_transform':[1,0,0],
      'template_metabolites':{},'template_reactions':{},'canvas':[width,height],
      'pathway_multiplicity':{'enabled':True,'normalization_flux':1},
      'circular_backbone':spec,
      'line_style':{'mode':'flux_width_only','threshold':.4,'width_range':[2.4,5.7],'width_exponent':.65,'unit':'relative flux (Pyruvate = 1)'},
      'publication_style':{'metabolite_font_pt':19,'reaction_font_pt':13,'legend_position':[width/2-210,height-100]},
      'hidden_metabolites':sorted(set(mids)&CURRENCY),'hidden_reactions':[],
      'metabolites':{mid:{'x':center[0],'y':center[1],'display_name':DISPLAY.get(mid,mid),'node_radius':5.5 if mid in CURRENCY else 9} for mid in mids},
      'reactions':{},
      'scientific_layout':{'enabled':True,'roles':{'source':['Formaldehyde'],'sink':['Pyruvate'],
        'cycle_anchor':['Ru5P'],'regeneration':['E4P','Xu5P','S7P','R5P','IML_TKT2','IML_TALA','IML_TKT1','IML_RPI','IML_RPE']},
        'role_palette':{'source':'#2F6F8F','sink':'#2F6F8F','cycle_anchor':'#8B6FA3','regeneration':'#B58AA6'},
        'reaction_styles':{rid:{'display_name':R_LABEL.get(rid,rid.removeprefix('IML_'))} for rid in reactions}}
    }
    for rid, reaction in reactions.items():
        cofactors = {}
        for mid in reaction['full_stoichiometry']:
            alias=ALIASES[mid]
            if alias in CURRENCY:
                cofactors[alias]={'x':center[0],'y':center[1],'display_name':DISPLAY.get(alias,alias),'node_radius':5.5,
                                 'scientific_role':'cofactor','color':'#9AA6B2'}
        layout['reactions'][rid]={'x':center[0],'y':center[1],
          'display_name':R_LABEL.get(rid,rid.removeprefix('IML_')),
          'pathway':'PPP' if rid in REGEN or rid=='IML_RPE' else 'ED' if rid in ('IML_EDD','IML_EDA') else 'Other',
          'color':'#B58AA6' if rid in REGEN or rid=='IML_RPE' else '#4A8D88' if rid not in ('RUMP_R05338','RUMP_R05339') else '#657786',
          'cofactors':cofactors}

    if index not in (1,4):
        cycle=index in (3,6,7,8,9)
        if index==2:
            seq=['DHAP','IML_MGSA','Methylglyoxal','IML_GLYOX3','D_lactate','IML_LDH_D','Pyruvate']
        elif index==5:
            seq=['DHAP','IML_TPI','G3P','IML_GAPD','BPG','IML_PGK','PG3','IML_PGM','PG2','IML_ENO','PEP','IML_PYK','Pyruvate']
        elif index in (3,6): seq=list(LOWER)
        elif index==7: seq=HP+['IML_HPYRRx']+TAIL
        elif index==8: seq=HP+['IML_HPYRRy']+TAIL
        else: seq=HP+['IML_HPYRI','TSA','IML_TRSARr']+TAIL
        n=len(seq)//2
        c2=[center[0]+1570,center[1]]
        r2=560
        if cycle:
            delta=360/n
            s2=180-delta*seq.index('IML_DHAPT')/2
        else:
            delta=250/n
            s2=150
        groups={}
        for i in range(n):
            source,rid=seq[2*i:2*i+2]
            target=seq[(2*i+2)%len(seq)]
            a0,a1=s2+i*delta,s2+(i+1)*delta
            amid=(a0+a1)/2
            a,b=point(c2,r2,a0),point(c2,r2,a1)
            side={'reaction':rid,'source':source,'target':target,'start':a,'end':b,
                  'control_points':bezier(c2,r2,a0,a1),'label_offset':[-68*math.cos(math.radians(amid)),-68*math.sin(math.radians(amid))-5]}
            for mid,point_,angle_,key in ((source,a,a0,'local_source'),(target,b,a1,'local_target')):
                local=mid in items[::2] or (index==6 and mid=='DHAP')
                if local:
                    side[key]=True
                    groups.setdefault(mid,{'metabolite':mid,'reactions':[],'position':point_,'label_offset':[0,0]})['reactions'].append(rid)
                    off=label_offset(DISPLAY.get(mid,mid),angle_)
                    groups[mid]['label_offset']=[off['label_dx'],off['label_dy']]
                    overrides[mid+'_'+rid]={'node_radius':9}
                else:
                    layout['metabolites'][mid].update(x=point_[0],y=point_[1])
                    overrides[mid]=label_offset(DISPLAY.get(mid,mid),angle_)
            spec['side_reactions'].append(side)
            # Position every nucleotide, redox cofactor and currency incidence.
            side_stoich={ALIASES[mid]:coef for mid,coef in reactions[rid]['full_stoichiometry'].items()}
            bysign={-1:[],1:[]}
            for mid,coef in side_stoich.items():
                if mid in CURRENCY: bysign[-1 if coef<0 else 1].append(mid)
            normal=[math.cos(math.radians(amid)),math.sin(math.radians(amid))]
            tangent=[-normal[1],normal[0]]
            for sign,metabolites in bysign.items():
                for j,mid in enumerate(metabolites):
                    along=sign*(65+j*85)
                    p=[c2[0]+(r2+140)*normal[0]+along*tangent[0],c2[1]+(r2+140)*normal[1]+along*tangent[1]]
                    cf=layout['reactions'][rid]['cofactors'][mid]
                    cf.update(x=p[0],y=p[1],align='center',**label_offset(DISPLAY.get(mid,mid),amid,radius=5.5))
        for mid,group in groups.items():
            group['reactions']=sorted(set(group['reactions']))
            if len(group['reactions'])>1:
                spec['shared_instances'][mid+'_output']=group
        if cycle:
            dha=[center[0]+860,center[1]+125]
            pyr=[center[0]+860,center[1]-150]
            layout['metabolites']['DHA'].update(x=dha[0],y=dha[1])
            layout['metabolites']['Pyruvate'].update(x=pyr[0],y=pyr[1])
            overrides['DHA']={'x':dha[0],'y':dha[1],'label_dx':0,'label_dy':48}
            overrides['Pyruvate']={'x':pyr[0],'y':pyr[1],'label_dx':0,'label_dy':-42}
            if index==6:
                dhap=[center[0]+650,center[1]+270]
                spec['shared_metabolites']['DHA']={'position':dha,'label_offset':[0,48]}
                spec['side_reactions'].append({'reaction':'IML_ALKP','source':'DHAP','target':'DHA',
                   'start':dhap,'end':dha,'label_offset':[10,38]})
                overrides['DHAP']={'x':dhap[0],'y':dhap[1],'label_dx':0,'label_dy':48}
                for j,mid in enumerate(('H2O','Pi')):
                    cf=layout['reactions']['IML_ALKP']['cofactors'][mid]
                    cf.update(x=dhap[0]+65+j*120,y=dhap[1]+160,align='center',label_dx=0,label_dy=38)
    if side_pathway_mode not in ('open_chain', 'legacy_circle'):
        raise ValueError('side_pathway_mode must be open_chain or legacy_circle')
    if index not in (1, 4) and side_pathway_mode == 'open_chain':
        from open_side_layout import unfold_side_chain
        side_ids = [s['reaction'] for s in spec['side_reactions']
                    if s['reaction'] not in ('IML_RPE', 'IML_ALKP')]
        stoich = {rid: {ALIASES[mid]: coefficient for mid, coefficient in r['full_stoichiometry'].items()}
                  for rid, r in reactions.items()}
        primary_mids = {mid for rid in items[1::2] for mid in stoich[rid] if mid not in CURRENCY}
        layout = unfold_side_chain(layout, side_ids, stoich, primary_mids)
    from publication_readability import enlarge_publication
    return enlarge_publication(layout)

def prepare():
    raw=SOURCE.read_bytes()
    assert hashlib.sha256(raw).hexdigest()==SOURCE_SHA, 'Archived source changed; do not silently mix runs.'
    pathways=[json.loads(line) for line in raw.decode('utf-8').splitlines() if line.strip()]
    assert len(pathways)==9
    (HERE/'source').mkdir(exist_ok=True)
    (HERE/'source/pathways.jsonl').write_bytes(raw)
    dump(HERE/'source/metabolite_aliases.json',ALIASES)
    from scripts_render_d3flux import read_model
    results=[]
    for i,p in enumerate(pathways,1):
        assert p['pathway_id']==f'pathway-{i:08d}'
        folder=HERE/f'candidate_{i:02d}'
        folder.mkdir(exist_ok=True)
        lines=[]
        for reaction in p['reactions']:
            sides=[]
            for sign in (-1,1):
                terms=[]
                for mid,coef in reaction['full_stoichiometry'].items():
                    if coef*sign>0:
                        terms.append(('' if abs(coef)==1 else f'{abs(coef):.12g}')+ALIASES[mid])
                sides.append(' + '.join(terms))
            lines.append('\t'.join([reaction['reaction_id'],*sides,'0']))
        (folder/'model.tsv').write_text('\n'.join(lines)+'\n',encoding='utf-8')
        dump(folder/'fit_result.json',{'opt_net_fluxes':{r['reaction_id']:r['flux'] for r in p['reactions']},
             'normalization':{'Pyruvate':1,'unit':'relative flux'},'pathway_id':p['pathway_id'],'source_sha256':SOURCE_SHA})
        dump(folder/'layout.json',make_layout(i,p))
        dump(folder/'source_pathway.json',p)
        model=read_model(folder/'model.tsv')
        for reaction in p['reactions']:
            actual={m.id:float(v) for m,v in model.reactions.get_by_id(reaction['reaction_id']).metabolites.items()}
            expected={ALIASES[mid]:float(v) for mid,v in reaction['full_stoichiometry'].items()}
            assert actual==expected,(i,reaction['reaction_id'])
        results.append({'candidate':i,'reactions':len(model.reactions),'reaction_and_stoichiometry_roundtrip':True,'flux_roundtrip':True})
    dump(HERE/'input_validation.json',{'source_sha256':SOURCE_SHA,'pathways':results})
    print(json.dumps(results,indent=2),flush=True)

def render(ids):
    import scripts_render_d3flux as native
    failed = []
    for i in ids:
        folder=HERE/f'candidate_{i:02d}'
        out=folder/f'candidate_{i:02d}'
        manifest_path=out.with_suffix('.manifest.json')
        previous_stamp=manifest_path.stat().st_mtime_ns if manifest_path.exists() else None
        cmd=[sys.executable,'-X','utf8',str(ROOT/'render_network.py'),str(folder/'fit_result.json'),str(folder/'model.tsv'),str(out.with_suffix('.html')),'--layout',str(folder/'layout.json')]
        sys.argv=cmd[3:]
        code=0
        with (folder/'render.log').open('w',encoding='utf-8') as log:
            with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
                try:
                    native.main()
                except Exception:
                    import traceback
                    traceback.print_exc()
                    code=1
        fresh_manifest=manifest_path.exists() and manifest_path.stat().st_mtime_ns != previous_stamp
        current_layout = json.loads((folder/'layout.json').read_text(encoding='utf-8'))
        if code and fresh_manifest and current_layout.get('publication_style', {}).get('auto_label_clearance'):
            report = json.loads(manifest_path.read_text(encoding='utf-8'))
            qa = report.get('rendered_elements', {}).get('strict_visual_qa', {})
            fatal = ('crossings', 'overlaps', 'node_hits', 'semantic_errors',
                     'path_continuity_errors', 'visual_connectivity_errors')
            label_errors = ('label_collisions', 'path_label_hits', 'clipped_labels', 'label_node_hits')
            if not any(qa.get(key) for key in fatal) and any(qa.get(key) for key in label_errors):
                from label_clearance import clear_publication_labels
                adjusted, placement = clear_publication_labels(current_layout, out.with_suffix('.html'))
                dump(folder/'label_clearance.json', placement)
                if placement['changes'] and not placement['unresolved']:
                    dump(folder/'layout.json', adjusted)
                    code = 0
                    with (folder/'render.log').open('a',encoding='utf-8') as log:
                        with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
                            print('\n=== Label-only clearance; repeat full native acceptance ===')
                            try:
                                native.main()
                            except Exception:
                                import traceback
                                traceback.print_exc()
                                code = 1
        print(f'Candidate {i:02d}: returncode={code}; svg={out.with_suffix(".svg").exists()}',flush=True)
        manifest=out.with_suffix('.manifest.json')
        if fresh_manifest:
            report=json.loads(manifest.read_text(encoding='utf-8'))
            qa=report['rendered_elements']['strict_visual_qa']
            print(json.dumps({k:v for k,v in qa.items() if k in ['pass','crossings','overlaps','label_collisions','path_label_hits','clipped_labels','node_hits','label_node_hits','semantic_errors']},ensure_ascii=False),flush=True)
        if code or not fresh_manifest or not report.get('acceptance', {}).get('pass'):
            failed.append(i)
    if failed:
        raise RuntimeError(f'Native render failed (stale manifests cannot pass): {failed}')

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('action',choices=['prepare','render'])
    parser.add_argument('--ids',default='1,2,3,4,5,6,7,8,9')
    args=parser.parse_args()
    prepare() if args.action=='prepare' else render([int(x) for x in args.ids.split(',')])

if __name__=='__main__': main()
