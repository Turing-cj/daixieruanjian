() => {
  // Audit rendered SVG geometry, not the planner's approximate paired graph.
  const svg=document.querySelector('#freeflux_map svg') || document.querySelector('svg');
  const step=4, paths=[...svg.querySelectorAll('path.link')];
  const distance=(a,b)=>Math.hypot(a.x-b.x,a.y-b.y);
  const sampled=paths.map(path=>{
    const length=path.getTotalLength(), n=Math.max(2,Math.ceil(length/step));
    const matrix=svg.getScreenCTM().inverse().multiply(path.getScreenCTM());
    return {path,d:path.__data__,points:Array.from({length:n+1},(_,i)=>{
      const p=path.getPointAtLength(length*i/n);
      return new DOMPoint(p.x,p.y).matrixTransform(matrix);
    })};
  });
  const name=x=>x.d.rxn.id+':'+(x.d.part || 'paired')+':'+x.d.source.id+'>'+x.d.target.id;
  const cross=(a,b,c)=>(b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x);
  function intersection(a,b,c,d) {
    if(cross(a,b,c)*cross(a,b,d)>=-1e-8 || cross(c,d,a)*cross(c,d,b)>=-1e-8) return null;
    const dx=b.x-a.x, dy=b.y-a.y, ex=d.x-c.x, ey=d.y-c.y;
    const t=((c.x-a.x)*ey-(c.y-a.y)*ex)/(dx*ey-dy*ex);
    return {x:a.x+t*dx,y:a.y+t*dy};
  }
  function pointSegment(p,a,b) {
    const dx=b.x-a.x,dy=b.y-a.y,l2=dx*dx+dy*dy;
    const t=l2 ? Math.max(0,Math.min(1,((p.x-a.x)*dx+(p.y-a.y)*dy)/l2)) : 0;
    return distance(p,{x:a.x+t*dx,y:a.y+t*dy});
  }
  function sharedJunctions(a,b) {
    const anchors=[];
    [a.d.source,a.d.target].forEach(node=>{
      if([b.d.source.id,b.d.target.id].includes(node.id)) anchors.push(node);
    });
    if(a.d.rxn.id===b.d.rxn.id) {
      const arc=a.d.rxn.notes.map_info.circle_arc;
      if(arc) [-1,1].forEach(sign=>{
        const theta=arc.reaction_angle_deg*Math.PI/180+sign*arc.join_half_length/arc.radius;
        anchors.push({x:arc.center[0]+arc.radius*Math.cos(theta),y:arc.center[1]+arc.radius*Math.sin(theta)});
      });
      const trunk=a.d.rxn.notes.map_info.local_trunk;
      if(trunk) {
        const local=a.d.part==='local' ? a : b.d.part==='local' ? b : null;
        const declared=[
          ...Object.values(trunk.input_joins || {}),
          ...Object.values(trunk.output_joins || {}),
        ];
        for(const p of declared) {
          // Exempt only a declared join close to the actual side trunk. Other
          // intersections of paths carrying the same reaction still fail.
          if(p && local && local.points.slice(1).some((q,i)=>pointSegment(p,local.points[i],q)<5)) anchors.push(p);
        }
      }
    }
    return anchors;
  }
  const crossings=[], overlaps=[];
  for(let i=0;i<sampled.length;i++) for(let j=i+1;j<sampled.length;j++) {
    const a=sampled[i],b=sampled[j], anchors=sharedJunctions(a,b);
    const allowed=p=>anchors.some(q=>distance(p,q)<24);
    let hit=false;
    for(let m=1;m<a.points.length && !hit;m++) for(let n=1;n<b.points.length;n++) {
      const p=intersection(a.points[m-1],a.points[m],b.points[n-1],b.points[n]);
      if(p && !allowed(p)) {hit=true;break;}
    }
    if(hit) crossings.push([name(a),name(b)]);
    let closeLength=0;
    for(const p of a.points) {
      if(allowed(p)) continue;
      if(b.points.slice(1).some((q,k)=>pointSegment(p,b.points[k],q)<1.25)) closeLength+=step;
    }
    if(closeLength>32) overlaps.push({paths:[name(a),name(b)],length_px:closeLength});
  }
  function visible(el) {
    const s=getComputedStyle(el);
    return s.display!=='none' && s.visibility!=='hidden' && +s.opacity!==0;
  }
  function box(el) {
    const b=el.getBBox(), matrix=svg.getScreenCTM().inverse().multiply(el.getScreenCTM());
    const a=new DOMPoint(b.x,b.y).matrixTransform(matrix);
    const c=new DOMPoint(b.x+b.width,b.y+b.height).matrixTransform(matrix);
    return {x:Math.min(a.x,c.x),y:Math.min(a.y,c.y),w:Math.abs(c.x-a.x),h:Math.abs(c.y-a.y)};
  }
  const labels=[...svg.querySelectorAll('text')].filter(visible).map(el=>({
    text:el.textContent,box:box(el)
  })).filter(l=>l.box.w>0 && l.box.h>0);
  const labelCollisions=[],labelHits=[],clipped=[];
  const vb=svg.viewBox.baseVal;
  for(let i=0;i<labels.length;i++) {
    const a=labels[i].box;
    if(a.x<0 || a.y<0 || a.x+a.w>vb.width || a.y+a.h>vb.height) clipped.push(labels[i].text);
    for(let j=i+1;j<labels.length;j++) {
      const b=labels[j].box;
      if(Math.min(a.x+a.w,b.x+b.w)-Math.max(a.x,b.x)>1 &&
         Math.min(a.y+a.h,b.y+b.h)-Math.max(a.y,b.y)>1) labelCollisions.push([labels[i].text,labels[j].text]);
    }
    for(const item of sampled) {
      if(item.points.some(p=>p.x>a.x+1 && p.x<a.x+a.w-1 && p.y>a.y+1 && p.y<a.y+a.h-1))
        labelHits.push([name(item),labels[i].text]);
    }
  }
  const nodes=[...svg.querySelectorAll('g.node-group')].map(g=>g.__data__).filter(Boolean);
  const labelNodeHits=[];
  for(const circle of [...svg.querySelectorAll('circle.metabolite')].filter(visible)) {
    const b=box(circle), center={x:b.x+b.w/2,y:b.y+b.h/2}, radius=Math.max(b.w,b.h)/2;
    for(const label of labels) {
      const a=label.box, nearest={x:Math.max(a.x,Math.min(center.x,a.x+a.w)),
                                 y:Math.max(a.y,Math.min(center.y,a.y+a.h))};
      if(distance(center,nearest)<radius+2) labelNodeHits.push([circle.id,label.text]);
    }
  }
  const nodeHits=[];
  for(const item of sampled) for(const node of nodes) {
    if([item.d.source.id,item.d.target.id,item.d.rxn.id].includes(node.id)) continue;
    if(item.points.some(p=>distance(p,node)<10)) nodeHits.push([name(item),node.id]);
  }
  const backbone=sampled.filter(x=>x.d.part==='backbone');
  let radialError=0;
  for(const item of backbone) {
    const arc=item.d.rxn.notes.map_info.circle_arc;
    for(const p of item.points) radialError=Math.max(radialError,Math.abs(Math.hypot(p.x-arc.center[0],p.y-arc.center[1])-arc.radius));
  }
  const snapshot=window.freefluxAuditSnapshot ? window.freefluxAuditSnapshot() : null;
  const semanticErrors=[], incidenceCounts={};
  const original=node=>node.notes.orig_id || node.id;
  for(const item of sampled) {
    const counts=incidenceCounts[item.d.rxn.id] || (incidenceCounts[item.d.rxn.id]={});
    const add=(node,sign)=>{const key=sign+':'+original(node);counts[key]=(counts[key]||0)+1;};
    if(['backbone','local','input'].includes(item.d.part)) add(item.d.source,-1);
    if(['backbone','local','output'].includes(item.d.part)) add(item.d.target,1);
    const marker=item.path.getAttribute('marker-end');
    const inputArrows=!!(snapshot && snapshot.notes.circular_backbone && snapshot.notes.circular_backbone.input_arrows);
    const needsArrow=['output','backbone','local'].includes(item.d.part) || item.d.part==='input' && inputArrows;
    if(item.d.part==='input' && !inputArrows && marker) semanticErrors.push('unexpected consumption arrow: '+name(item));
    if(needsArrow && !marker) semanticErrors.push('missing directed incidence arrow: '+name(item));
    if(marker) {
      const match=marker.match(/^url\(#(.+)\)$/), element=match && document.getElementById(match[1]);
      if(!element || element.tagName.toLowerCase()!=='marker' || !match[1].endsWith(item.d.rxn.id))
        semanticErrors.push('invalid or wrong-reaction arrow marker: '+name(item));
    }
  }
  if(snapshot && snapshot.notes.circular_backbone) for(const reaction of snapshot.reactions) {
    const info=reaction.notes.map_info;
    if(info.hidden) continue;
    const reverse=info.flux < -1E-10 ? -1 : 1, expected={};
    for(const [id,c] of Object.entries(reaction.metabolites)) expected[(Math.sign(c)*reverse)+':'+id]=1;
    const got=incidenceCounts[reaction.id] || {};
    if(JSON.stringify(Object.entries(expected).sort())!==JSON.stringify(Object.entries(got).sort()))
      semanticErrors.push({reaction:reaction.id,expected,got});
  }
  const presentation=snapshot && snapshot.notes.circular_backbone;
  const nodeGroups=[...svg.querySelectorAll('g.node-group')].filter(g=>g.__data__ && g.__data__.type!=='rxn');
  const pathContinuityErrors=[], visualConnectivityErrors=[];
  const visibleNodes=new Map(nodeGroups.map(g=>[g.__data__.id,g]));
  // A complete incidence table is insufficient: it used to accept a path
  // containing multiple disconnected SVG subpaths, or a detached display copy.
  for(const item of sampled) {
    const d=item.path.getAttribute('d') || '', style=getComputedStyle(item.path);
    if((d.match(/[mM]/g) || []).length!==1 || !Number.isFinite(item.path.getTotalLength()) ||
       item.path.getTotalLength()<=0 || !visible(item.path) || style.stroke==='none' || +style.strokeOpacity===0)
      pathContinuityErrors.push({path:name(item),reason:'missing, hidden or fragmented stroke'});
    const stroke=parseFloat(style.strokeWidth) || 1;
    for(const [node,point,hasArrow] of [
      [['backbone','local','input'].includes(item.d.part) ? item.d.source : null,item.points[0],false],
      [['backbone','local','output'].includes(item.d.part) ? item.d.target : null,item.points[item.points.length-1],true]]) {
      if(!node || node.type==='rxn') continue;
      const group=visibleNodes.get(node.id);
      if(!group || group.__data__!==node) {
        visualConnectivityErrors.push({path:name(item),reason:'endpoint has no matching rendered node',node:node.id});
        continue;
      }
      const circle=group.querySelector('circle.metabolite');
      if(!circle || !visible(circle)) {
        pathContinuityErrors.push({path:name(item),reason:'endpoint node is not visible',node:node.id});
        continue;
      }
      const b=box(circle), center={x:b.x+b.w/2,y:b.y+b.h/2}, radius=Math.max(b.w,b.h)/2;
      // Include the deliberate node clearance and arrowhead; never permit a
      // long disconnected end merely because its biochemical ID is correct.
      const allowance=radius+7+(hasArrow ? Math.max(14,stroke*3) : 0);
      if(distance(point,center)>allowance)
        pathContinuityErrors.push({path:name(item),reason:'stroke endpoint detached from node',node:node.id,
                                   gap_px:distance(point,center)-radius});
    }
    if(item.d.part==='input' || item.d.part==='output') {
      const trunk=sampled.find(p=>p.d.rxn.id===item.d.rxn.id && ['backbone','local'].includes(p.d.part));
      const point=item.d.part==='input' ? item.points[item.points.length-1] : item.points[0];
      const allowance=item.d.part==='input' ? Math.max(18,stroke*3) : 4;
      if(!trunk || !trunk.points.slice(1).some((q,i)=>pointSegment(point,trunk.points[i],q)<=allowance))
        pathContinuityErrors.push({path:name(item),reason:'reaction branch detached from its trunk'});
    }
    if(item.points.some(p=>p.x<0 || p.y<0 || p.x>vb.width || p.y>vb.height))
      pathContinuityErrors.push({path:name(item),reason:'stroke outside canvas'});
  }
  function unionFind() {
    const parent=new Map();
    const root=x=>{if(!parent.has(x)) parent.set(x,x); if(parent.get(x)!==x) parent.set(x,root(parent.get(x))); return parent.get(x);};
    return {root,join:(a,b)=>parent.set(root(a),root(b)),keys:()=>[...parent.keys()]};
  }
  const chemical=unionFind(), visual=unionFind();
  const carbon=node=>node.type!=='rxn' && (!node.cofactor || node.local_instance);
  if(presentation && snapshot) {
    const carbonIds=new Set(snapshot.metabolites.filter(m=>!m.notes.map_info.hidden).map(m=>m.id));
    for(const r of snapshot.reactions) {
      if(r.notes.map_info.hidden) continue;
      const mids=Object.keys(r.metabolites).filter(mid=>carbonIds.has(mid));
      if(!mids.length) continue;
      visual.root('r:'+r.id);
      for(const mid of mids) chemical.join('r:'+r.id,'m:'+mid);
    }
    for(const item of sampled) for(const node of [item.d.source,item.d.target])
      if(carbon(node)) visual.join('r:'+item.d.rxn.id,'v:'+node.id);
    for(const g of nodeGroups) if(carbon(g.__data__)) visual.root('v:'+g.__data__.id);
    const chemicalGroups=new Map();
    for(const key of visual.keys()) {
      const g=key.startsWith('v:') ? visibleNodes.get(key.slice(2)) : null;
      const logical=key.startsWith('r:') ? key : g ? 'm:'+original(g.__data__) : null;
      if(!logical) continue;
      const component=chemical.root(logical);
      if(!chemicalGroups.has(component)) chemicalGroups.set(component,new Map());
      const groups=chemicalGroups.get(component), connected=visual.root(key);
      if(!groups.has(connected)) groups.set(connected,[]);
      groups.get(connected).push(key);
    }
    for(const groups of chemicalGroups.values()) if(groups.size>1)
      visualConnectivityErrors.push({reason:'one biochemical carbon component is split into disconnected visual components',
                                     components:[...groups.values()].map(ids=>ids.sort())});
  }
  const multiplicity=snapshot && snapshot.notes.pathway_multiplicity;
  const localAmounts=!!(multiplicity && multiplicity.label_mode==='local_flow');
  const sourceSplit=!!(multiplicity && multiplicity.relationship_mode==='source_split');
  if(presentation) {
    // Original chemical IDs alone cannot prove visual connectivity: an
    // unconnected clone has the same original ID. Check actual node identity.
    for(const mid of presentation.shared_metabolites || []) {
      const views=nodeGroups.filter(g=>original(g.__data__)===mid);
      if(views.length!==1 || views[0].__data__.id!==mid) semanticErrors.push('shared pool is split: '+mid);
      const canonical=views.length===1 ? views[0].__data__ : null;
      for(const item of sampled) for(const endpoint of [item.d.source,item.d.target])
        if(original(endpoint)===mid && endpoint!==canonical) semanticErrors.push('disconnected shared pool incidence: '+name(item));
    }
    for(const [viewId,group] of Object.entries(presentation.shared_instances || {})) {
      const views=nodeGroups.filter(g=>g.__data__.id===viewId);
      if(views.length!==1 || original(views[0].__data__)!==group.metabolite) semanticErrors.push('invalid shared instance: '+viewId);
      const canonical=views.length===1 ? views[0].__data__ : null;
      for(const rid of group.reactions) {
        const ends=sampled.filter(p=>p.d.rxn.id===rid).flatMap(p=>[p.d.source,p.d.target]).filter(n=>original(n)===group.metabolite);
        if(ends.length!==1 || ends[0]!==canonical) semanticErrors.push('shared instance incidence is split: '+viewId+'/'+rid);
      }
    }
    for(const [mid,code] of Object.entries(presentation.pool_codes || {})) {
      if(sourceSplit) continue; // Independently check source/split labels below.
      const info=snapshot.metabolites.find(m=>m.id===mid).notes.map_info;
      const expected=(info.instance_base_name || mid)+' ['+code+']';
      for(const g of nodeGroups.filter(g=>original(g.__data__)===mid)) {
        const label=[...g.querySelectorAll('text.nodelabel tspan.text')].map(t=>t.textContent).join('\n');
        const base=info.pathway_multiplicity===undefined ? label : label.replace(/^\S+ /,'');
        if(base!==expected) semanticErrors.push('incorrect pool label: '+g.__data__.id+' / '+label);
      }
    }
  }
  if(multiplicity && multiplicity.enabled) {
    // Independent recomputation from chemistry: agreeing copies of a wrong
    // cached value must not pass a publication audit.
    const scale=multiplicity.normalization_flux, totals={};
    if(!Number.isFinite(scale) || scale<=0) semanticErrors.push('invalid multiplicity normalization');
    for(const reaction of snapshot.reactions) {
      const flux=reaction.notes.map_info.flux;
      if(!Number.isFinite(flux)) semanticErrors.push('invalid multiplicity flux: '+reaction.id);
      for(const [mid,coefficient] of Object.entries(reaction.metabolites)) {
        const value=coefficient*flux/scale, pair=totals[mid] || (totals[mid]=[0,0]);
        pair[value>=0 ? 0 : 1]+=Math.abs(value);
      }
    }
    const localTotals={}, localIncidences={}, seenIncidences=new Set();
    if(localAmounts) {
      const reactions=new Map(snapshot.reactions.map(r=>[r.id,r]));
      for(const item of sampled) {
        const ends=[];
        if(['backbone','local','input'].includes(item.d.part)) ends.push(item.d.source);
        if(['backbone','local','output'].includes(item.d.part)) ends.push(item.d.target);
        for(const node of ends) {
          const mid=original(node), reaction=reactions.get(item.d.rxn.id);
          const key=JSON.stringify([item.d.rxn.id,mid]);
          if(seenIncidences.has(key)) continue; // duplicates independently fail incidenceCounts above
          seenIncidences.add(key);
          const value=reaction.metabolites[mid]*reaction.notes.map_info.flux/scale;
          const pair=localTotals[node.id] || (localTotals[node.id]=[0,0]);
          pair[value>=0 ? 0 : 1]+=Math.abs(value);
          (localIncidences[node.id] || (localIncidences[node.id]=[])).push({reaction:reaction.id,signed_rate:value});
        }
      }
      for(const [mid,pair] of Object.entries(totals)) {
        const expected=Math.max(...pair), cached=multiplicity.counts && multiplicity.counts[mid];
        if(!Number.isFinite(cached) || Math.abs(cached-expected)>1e-9*Math.max(1,expected))
          semanticErrors.push('incorrect separate pool total: '+mid);
      }
    }
    for(const g of nodeGroups) {
      const mid=original(g.__data__), info=g.__data__.notes.map_info;
      const poolTotal=Math.max(...(totals[mid] || [0,0]));
      const pair=(localAmounts ? localTotals[g.__data__.id] : totals[mid]) || [0,0];
      const expected=Math.max(...pair);
      const label=[...g.querySelectorAll('text.nodelabel tspan.text')].map(t=>t.textContent).join('\n');
      const match=label.match(/^(\S+) (.+)$/), shown=match ? Number(match[1]) : NaN;
      const tolerance=1e-9*Math.max(1,expected);
      if(!Number.isFinite(expected) || !Number.isFinite(shown) || Math.abs(shown-expected)>tolerance ||
         !Number.isFinite(info.pathway_multiplicity) || Math.abs(info.pathway_multiplicity-expected)>tolerance ||
         !match || match[2]!==info.pathway_label_base)
        semanticErrors.push('incorrect pathway multiplicity: '+g.__data__.id+' / '+label);
      if(localAmounts && (info.amount_scope!=='local_flow' ||
          !Number.isFinite(info.pool_total) || Math.abs(info.pool_total-poolTotal)>1e-9*Math.max(1,poolTotal)))
        semanticErrors.push('incorrect local-flow scope or pool total: '+g.__data__.id);
    }
    if(sourceSplit) {
      // This recomputation uses actual link endpoints and signed chemistry,
      // never the Python provenance cache or a node's claimed parent.
      const same=(a,b)=>JSON.stringify(a)===JSON.stringify(b);
      const near=(a,b)=>Number.isFinite(a) && Math.abs(a-b)<=1e-9*Math.max(1,Math.abs(b));
      const compare=(a,b)=>a<b ? -1 : a>b ? 1 : 0;
      const ledger=multiplicity.source_split;
      if(!localAmounts || !ledger || ledger.version!==1)
        semanticErrors.push('invalid source/split contract');
      for(const [mid,code] of Object.entries(presentation && presentation.pool_codes || {})) {
        const groups=nodeGroups.filter(g=>original(g.__data__)===mid)
          .sort((a,b)=>compare(a.__data__.id,b.__data__.id));
        const sources=groups.flatMap(g=>(localIncidences[g.__data__.id] || [])
          .filter(i=>i.signed_rate>0).map(i=>({reaction:i.reaction,view_id:g.__data__.id,amount:i.signed_rate})))
          .sort((a,b)=>compare(a.reaction,b.reaction)||compare(a.view_id,b.view_id));
        sources.forEach((s,i)=>{s.id=sources.length===1 ? code : code+(i+1);});
        const ids=sources.map(s=>s.id), total=totals[mid] || [0,0];
        const tolerance=1e-10*Math.max(1,...total);
        const partition=sources.length===1 && total[1]<=total[0]+tolerance;
        const children=groups.filter(g=>{
          const pair=localTotals[g.__data__.id] || [0,0]; return pair[0]===0 && pair[1]>0;
        }).map(g=>g.__data__.id);
        const uses=[];
        for(const g of groups) {
          const vid=g.__data__.id, info=g.__data__.notes.map_info;
          const pair=localTotals[vid] || [0,0], own=sources.filter(s=>s.view_id===vid);
          let role,marker,parents=[];
          const poolMarker=ids.length>1 ? ids.join('+')+':pool' : code+':pool';
          if(own.length && pair[1]>pair[0]+tolerance) {role='source_and_pool_use';marker=poolMarker;}
          else if(own.length) {role='source';marker=own.map(s=>s.id).join('+');}
          else if(partition && children.includes(vid)) {role='split_use';marker=ids[0]+'.'+(children.indexOf(vid)+1);parents=ids;}
          else if(pair[1]>0) {role='pool_use';marker=poolMarker;}
          else {role='inactive';marker=code+':0';}
          const expected={metabolite:mid,pool_code:code,role,marker,
            source_ids:own.map(s=>s.id),source_reactions:own.map(s=>s.reaction),
            parent_source_ids:parents,possible_source_ids:ids};
          const originalInfo=snapshot.metabolites.find(m=>m.id===mid).notes.map_info;
          const expectedBase=(originalInfo.instance_base_name || mid)+' ['+marker+']';
          const label=[...g.querySelectorAll('text.nodelabel tspan.text')].map(t=>t.textContent).join('\n');
          if(label.replace(/^\S+ /,'')!==expectedBase)
            semanticErrors.push('incorrect source/split label: '+vid+' / '+label);
          for(const actual of [info.flow_relation, ledger && ledger.views && ledger.views[vid]])
            if(!actual || Object.keys(expected).some(key=>!same(actual[key],expected[key])))
              semanticErrors.push('incorrect source/split metadata: '+vid);
          for(const inc of [...(localIncidences[vid] || [])].sort((a,b)=>compare(a.reaction,b.reaction)))
            if(inc.signed_rate<0) uses.push({view_id:vid,reaction:inc.reaction,amount:-inc.signed_rate,
              parent_source_ids:partition ? ids : [],possible_source_ids:ids});
        }
        const pool=ledger && ledger.pools && ledger.pools[mid];
        const recordsMatch=(actual,expected,fields)=>Array.isArray(actual) && actual.length===expected.length &&
          expected.every((e,i)=>fields.every(key=>same(actual[i][key],e[key])) && near(actual[i].amount,e.amount));
        if(!pool || pool.single_source_partition!==partition ||
           pool.boundary_import_required!==(total[1]>total[0]+tolerance) ||
           !near(pool.produced,total[0]) || !near(pool.consumed,total[1]) ||
           !recordsMatch(pool.sources,sources,['id','reaction','view_id']) ||
           !recordsMatch(pool.consumption_allocations,uses,['view_id','reaction','parent_source_ids','possible_source_ids']))
          semanticErrors.push('incorrect source/split source-allocation ledger: '+mid);
      }
      const legend='[A.1], [A.2] = portions of [A]; [B1], [B2] = independent sources; [B1+B2:pool] = pooled use from those possible sources; exact allocation unresolved.';
      if(![...svg.querySelectorAll('.circle-caption text')].some(t=>t.textContent===legend))
        semanticErrors.push('missing source/split legend');
    }
  }
  const circular=!!(snapshot && snapshot.notes.circular_backbone);
  const expectedArcs=circular ? snapshot.notes.circular_backbone.backbone_reaction_count : null;
  const screenMatrix=svg.getScreenCTM();
  const screenAxisRatio=Math.hypot(screenMatrix.a,screenMatrix.b)/Math.hypot(screenMatrix.c,screenMatrix.d);
  const pass=circular && backbone.length===expectedArcs && radialError<0.05 && Math.abs(screenAxisRatio-1)<0.001 &&
    !crossings.length && !overlaps.length && !labelCollisions.length && !labelHits.length &&
    !clipped.length && !nodeHits.length && !labelNodeHits.length && !semanticErrors.length &&
    !pathContinuityErrors.length && !visualConnectivityErrors.length;
  return {pass,sample_step_px:step,backbone_arc_count:backbone.length,
    quantity_label_mode:localAmounts ? 'local_flow' : 'pool_total',
    relationship_mode:sourceSplit ? 'source_split' : 'pool',
    screen_axis_ratio:screenAxisRatio,
    maximum_rendered_radial_error_px:radialError,crossings,overlaps,
    label_collisions:labelCollisions,path_label_hits:labelHits,clipped_labels:clipped,
    node_hits:nodeHits,label_node_hits:labelNodeHits,semantic_errors:semanticErrors,
    path_continuity_errors:pathContinuityErrors,visual_connectivity_errors:visualConnectivityErrors,
    visual_carbon_component_count:new Set(visual.keys().map(visual.root)).size,
    biochemical_carbon_component_count:new Set(chemical.keys().map(chemical.root)).size,
    explicit_input_arrow_count:sampled.filter(p=>p.d.part==='input' && p.path.getAttribute('marker-end')).length,
    shared_pool_count:(presentation && presentation.shared_metabolites || []).length,
    shared_instance_count:Object.keys(presentation && presentation.shared_instances || {}).length,
    rendered_reaction_count:Object.keys(incidenceCounts).length,
    snapshot_available:!!snapshot};
}
