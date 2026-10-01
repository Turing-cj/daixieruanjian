() => {
  // Move labels only, using the final browser font metrics. No node, reaction,
  // path or arrow may be removed, shortened, or moved to create whitespace.
  const svg=document.querySelector('#freeflux_map svg'), vb=svg.viewBox.baseVal;
  const shown=e=>{const s=getComputedStyle(e);return s.display!=='none' && s.visibility!=='hidden' && +s.opacity!==0;};
  function box(e) {
    const b=e.getBBox(), m=svg.getScreenCTM().inverse().multiply(e.getScreenCTM());
    const p=new DOMPoint(b.x,b.y).matrixTransform(m), q=new DOMPoint(b.x+b.width,b.y+b.height).matrixTransform(m);
    return {x:Math.min(p.x,q.x),y:Math.min(p.y,q.y),w:Math.abs(q.x-p.x),h:Math.abs(q.y-p.y)};
  }
  const paths=[...svg.querySelectorAll('path.link')], before=paths.map(p=>p.getAttribute('d'));
  const geometry=paths.map(p=>{
    const n=Math.max(2,Math.ceil(p.getTotalLength()/3));
    const m=svg.getScreenCTM().inverse().multiply(p.getScreenCTM());
    return {box:box(p),pad:Math.max(3,(parseFloat(getComputedStyle(p).strokeWidth)||1)/2+2),
      points:Array.from({length:n+1},(_,i)=>{const q=p.getPointAtLength(p.getTotalLength()*i/n);return new DOMPoint(q.x,q.y).matrixTransform(m);})};
  });
  const circles=[...svg.querySelectorAll('circle.metabolite')].filter(shown).map(e=>{
    const b=box(e);return {x:b.x+b.w/2,y:b.y+b.h/2,r:Math.max(b.w,b.h)/2+4};
  });
  const labels=[...svg.querySelectorAll('text')].filter(shown).map(e=>({e,b:box(e),
    id:e.getAttribute('data-label-id'),data:e.__data__})).filter(l=>l.b.w>0 && l.b.h>0);
  const intersect=(a,b,p=0)=>a.x+a.w+p>b.x && b.x+b.w+p>a.x && a.y+a.h+p>b.y && b.y+b.h+p>a.y;
  function clear(rect,self) {
    if(rect.x<4 || rect.y<4 || rect.x+rect.w>vb.width-4 || rect.y+rect.h>vb.height-4) return false;
    for(const other of labels) if(other!==self && intersect(rect,other.b,4)) return false;
    for(const c of circles) {
      const x=Math.max(rect.x,Math.min(c.x,rect.x+rect.w)), y=Math.max(rect.y,Math.min(c.y,rect.y+rect.h));
      if(Math.hypot(x-c.x,y-c.y)<c.r) return false;
    }
    for(const path of geometry) {
      if(!intersect(rect,path.box,path.pad)) continue;
      if(path.points.some(p=>p.x>rect.x-path.pad && p.x<rect.x+rect.w+path.pad &&
                             p.y>rect.y-path.pad && p.y<rect.y+rect.h+path.pad)) return false;
    }
    return true;
  }
  const shifts=[[0,0]];
  for(let radius=8;radius<=192;radius+=8) for(let i=0;i<24;i++)
    shifts.push([radius*Math.cos(2*Math.PI*i/24),radius*Math.sin(2*Math.PI*i/24)]);
  const changes=new Map(), movable=labels.filter(l=>l.id && l.data && l.data.notes && l.data.notes.map_info)
    .sort((a,b)=>a.id.localeCompare(b.id));
  for(let pass=0;pass<4;pass++) {
    let moved=0;
    for(const label of movable) {
      if(clear(label.b,label)) continue;
      const origin=label.b;
      for(const [dx,dy] of shifts.slice(1)) {
        const candidate={...origin,x:origin.x+dx,y:origin.y+dy};
        if(!clear(candidate,label)) continue;
        const info=label.data.notes.map_info;
        info.label_dx=(+info.label_dx||0)+dx; info.label_dy=(+info.label_dy||0)+dy;
        label.e.setAttribute('transform','translate('+info.label_dx+','+info.label_dy+')');
        label.b=box(label.e);
        changes.set(label.id,{id:label.id,label_dx:info.label_dx,label_dy:info.label_dy,text:label.e.textContent});
        moved++; break;
      }
    }
    if(!moved) break;
  }
  return {changes:[...changes.values()], unresolved:movable.filter(l=>!clear(l.b,l)).map(l=>l.id),
          paths_unchanged:paths.every((p,i)=>p.getAttribute('d')===before[i])};
}
