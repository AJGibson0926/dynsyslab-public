// Self-contained interactive HTML of a Viewer payload (PlotSvg.hpp payloadHtml). The page holds the
// plotted samples (already in the figure's value scale and projection) and a short script that
// draws them on HTML canvases (WHATWG HTML Living Standard, canvas 2D context); it needs no network.
#include "PlotSvg.hpp"

#include "dynsyslab_viewer/PlotScale.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>

namespace dsl {
namespace {
namespace json = dynsyslab::core::json;

json::Value numbers(const std::vector<double>& values) {
    json::Value::Array out;
    out.reserve(values.size());
    for (const double value : values) out.emplace_back(std::isfinite(value) ? json::Value(value) : json::Value{});
    return json::Value(std::move(out));
}

json::Value traceJson(const Figure& figure, const Trace& trace) {
    json::Value::Object out{{"label", json::Value(plainLabel(trace.label))}, {"color", json::Value(trace.color)},
                            {"style", json::Value(trace.lineStyle)},
                            {"marker", json::Value(trace.marker || trace.lineStyle == "markers")}};
    if (trace.rightAxis) out["right"] = json::Value(true);  // reads off the figure's right ordinate
    if (trace.heatmap && trace.zShape.size() == 2) {
        out["heat"] = numbers(trace.z);
        out["rows"] = json::Value(static_cast<double>(trace.zShape[0]));
        out["cols"] = json::Value(static_cast<double>(trace.zShape[1]));
        return json::Value(std::move(out));
    }
    const auto p = projectedTrace(figure, trace);
    std::vector<double> y = p.y;
    if (figure.kind != "three_d") {
        const auto scale = ordinateScale(figure);
        for (auto& value : y) value = scale.forward(value);
    }
    out["x"] = numbers(p.x);
    out["y"] = numbers(y);
    return json::Value(std::move(out));
}

json::Value payloadJson(const Payload& payload) {
    json::Value::Array pages;
    for (int page = 0; page < pageCount(payload); ++page) {
        json::Value::Array figures;
        for (const auto index : pageFigures(payload, page)) {
            const auto& figure = payload.figures[index];
            json::Value::Array traces;
            for (const auto& trace : figure.traces) traces.push_back(traceJson(figure, trace));
            json::Value::Object item{{"title", json::Value(plainLabel(figure.title))}, {"x_label", json::Value(plainLabel(figure.xLabel))},
                {"y_label", json::Value(plainLabel(figure.yLabel) + (figure.kind == "three_d" ? "" : ordinateScale(figure).kind
                    == ValueScaleKind::logarithmic ? " (log10)" : ordinateScale(figure).kind
                    == ValueScaleKind::symmetric_log ? " (symlog)" : ""))},
                {"traces", json::Value(std::move(traces))}};
            if (figure.focusX) item["focus"] = numbers({figure.focusX->first, figure.focusX->second});
            item["kind"] = json::Value(figure.kind);
            item["legend"] = json::Value(figure.legendPlacement);
            // The right ordinate (value = (plotted - offset) / scale) and stacked bands, on linear plots.
            if (figure.rightAxis && figure.kind != "three_d" && ordinateScale(figure).kind == ValueScaleKind::linear) {
                const auto& axis = *figure.rightAxis;
                json::Value::Object right{{"label", json::Value(plainLabel(axis.label))},
                                          {"scale", json::Value(axis.scale)}, {"offset", json::Value(axis.offset)}};
                if (axis.leftBand && axis.rightBand) {
                    right["left_band"] = numbers({axis.leftBand->first, axis.leftBand->second});
                    right["right_band"] = numbers({axis.rightBand->first, axis.rightBand->second});
                }
                item["right"] = json::Value(std::move(right));
            }
            figures.emplace_back(std::move(item));
        }
        const auto label = static_cast<std::size_t>(page) < payload.pageLabels.size()
            ? payload.pageLabels[static_cast<std::size_t>(page)] : "Page " + std::to_string(page + 1);
        pages.emplace_back(json::Value::Object{{"label", json::Value(label)}, {"figures", json::Value(std::move(figures))}});
    }
    return json::Value(json::Value::Object{{"columns", json::Value(static_cast<double>(std::max(1, payload.columns)))},
                                           {"pages", json::Value(std::move(pages))}});
}

constexpr const char* kScript = R"JS(
const D=JSON.parse(document.getElementById('dsl-data').textContent),P=['#1F77B4','#D62728','#2CA02C','#FF7F0E','#9467BD','#8C564B','#E377C2','#7F7F7F'];
const root=document.getElementById('dsl-root'),tabs=document.getElementById('dsl-tabs');
function nice(a,b){const r=(b-a)/5,m=Math.pow(10,Math.floor(Math.log10(r))),f=r/m,s=m*(f<1.5?1:f<3.5?2:f<7.5?5:10),t=[];for(let v=Math.ceil(a/s)*s;v<=b+1e-9*(b-a);v+=s)t.push(Math.abs(v)<1e-12*(b-a)?0:v);return t;}
function fmt(v){return Math.abs(v)>=1e4||(Math.abs(v)<1e-3&&v!==0)?v.toExponential(2):(+v.toPrecision(4)).toString();}
function viridis(t){const s=[[68,1,84],[59,82,139],[33,145,140],[94,201,98],[253,231,37]];t=Math.min(1,Math.max(0,t||0))*4;const k=Math.min(3,Math.floor(t)),f=t-k;return 'rgb('+s[k].map((c,i)=>Math.round(c+f*(s[k+1][i]-c))).join(',')+')';}
function limits(f){let x0=1/0,x1=-1/0,y0=1/0,y1=-1/0;for(const t of f.traces){if(!t.x)continue;for(let i=0;i<t.x.length;i++){const x=t.x[i],y=t.y[i];if(x==null||y==null)continue;if(f.focus&&(x<f.focus[0]||x>f.focus[1]))continue;x0=Math.min(x0,x);x1=Math.max(x1,x);y0=Math.min(y0,y);y1=Math.max(y1,y);}}
 if(f.focus){x0=f.focus[0];x1=f.focus[1];}if(!isFinite(x0)){x0=0;x1=1;y0=0;y1=1;}if(x1<=x0){x0-=.5;x1+=.5;}if(y1<=y0){y0-=.5;y1+=.5;}const p=.05*(y1-y0);return [x0,x1,y0-p,y1+p];}
function panel(f){const box=document.createElement('div');box.className='fig';const c=document.createElement('canvas'),o=document.createElement('div');o.className='readout';box.append(c,o);
 let L=limits(f);const home=L.slice();const R=f.right,W=()=>c.clientWidth,H=()=>c.clientHeight,m={l:64,r:R?78:16,t:30,b:44};
 const sx=x=>m.l+(x-L[0])/(L[1]-L[0])*(W()-m.l-m.r),sy=y=>H()-m.b-(y-L[2])/(L[3]-L[2])*(H()-m.t-m.b);
 const ix=p=>L[0]+(p-m.l)/(W()-m.l-m.r)*(L[1]-L[0]),iy=p=>L[2]+(H()-m.b-p)/(H()-m.t-m.b)*(L[3]-L[2]);
 const entries=f.traces.map((t,k)=>({t,col:t.color||P[k%P.length]})).filter(e=>e.t.label&&!e.t.heat).slice(0,8);
 let place=null;
 function dash(g,s){g.setLineDash(s==='dash'||s==='time-dash'?[7,4]:s==='dot'?[1.5,3]:s==='dash-dot'?[7,3,1.5,3]:[]);}
 function sample(g,e,x,y){g.strokeStyle=g.fillStyle=e.col;if(e.t.marker){g.setLineDash([]);g.beginPath();g.arc(x+9,y,3,0,7);g.fill();}else{g.lineWidth=2;dash(g,e.t.style);g.beginPath();g.moveTo(x,y);g.lineTo(x+18,y);g.stroke();}}
 // Legend placement, decided once on the home view: inside, in the first corner with no data
 // under it; a time or spectrum plot with data in every corner, or "top", gets one row under the title.
 function decide(g){if(!entries.length)return {kind:'none'};g.font='11px Segoe UI,Helvetica,sans-serif';
  const bw=Math.max(...entries.map(e=>26+g.measureText(e.t.label).width)),bh=16*entries.length+4;if(f.legend==='top')return {kind:'row'};
  const pw=W()-m.l-m.r,ph=H()-m.t-m.b,cs=[[m.l+pw-bw-6,m.t+6],[m.l+6,m.t+6],[m.l+pw-bw-6,m.t+ph-bh-6],[m.l+6,m.t+ph-bh-6]],searched=f.kind==='time'||f.kind==='fourier';
  for(const q of cs){let clear=true;if(searched)for(const t of f.traces){if(!t.x||!clear)continue;for(let i=0;i<t.x.length;i++){if(t.x[i]==null||t.y[i]==null)continue;const px=sx(t.x[i]),py=sy(t.y[i]);if(px>=q[0]-4&&px<=q[0]+bw+4&&py>=q[1]-4&&py<=q[1]+bh+4){clear=false;break;}}}
   if(clear)return {kind:'corner',x:q[0],y:q[1]};}
  return f.legend?{kind:'corner',x:cs[0][0],y:cs[0][1]}:{kind:'row'};}
 function legend(g){g.font='11px Segoe UI,Helvetica,sans-serif';g.textAlign='left';
  if(place.kind==='row'){const ws=entries.map(e=>26+g.measureText(e.t.label).width+18);let x=Math.max(8,(W()-ws.reduce((a,b)=>a+b,0)+18)/2);
   entries.forEach((e,k)=>{sample(g,e,x,38);g.fillStyle='#111';g.fillText(e.t.label,x+24,42);x+=ws[k];});}
  else if(place.kind==='corner')entries.forEach((e,k)=>{const y=place.y+10+16*k;sample(g,e,place.x+2,y);g.fillStyle='#111';g.fillText(e.t.label,place.x+26,y+4);});
  g.setLineDash([]);}
 function yTick(g,at,text,right){g.beginPath();g.moveTo(m.l,sy(at));g.lineTo(W()-m.r,sy(at));g.stroke();g.textAlign=right?'left':'right';g.fillText(text,right?W()-m.r+5:m.l-5,sy(at)+4);}
 function yTitle(g,text,lo,hi,right){if(!text)return;const a=Math.max(L[2],Math.min(lo,hi)),b=Math.min(L[3],Math.max(lo,hi));if(!(b>a))return;g.save();g.translate(right?W()-m.r+62:14,(sy(a)+sy(b))/2);g.rotate(-Math.PI/2);g.textAlign='center';g.fillText(text,0,0);g.restore();}
 function draw(){const r=window.devicePixelRatio||1;c.width=W()*r;c.height=H()*r;const g=c.getContext('2d');g.setTransform(r,0,0,r,0,0);g.clearRect(0,0,W(),H());
  if(!place){const saved=L;L=home;place=decide(g);L=saved;if(place.kind==='row')m.t=50;}
  g.font='13px Segoe UI,Helvetica,sans-serif';g.fillStyle='#111';g.textAlign='center';g.fillText(f.title,W()/2,18);g.font='11px Segoe UI,Helvetica,sans-serif';
  g.strokeStyle='#e5e5e5';for(const t of nice(L[0],L[1])){g.beginPath();g.moveTo(sx(t),m.t);g.lineTo(sx(t),H()-m.b);g.stroke();g.fillText(fmt(t),sx(t),H()-m.b+14);}
  if(R&&R.left_band&&R.right_band){  // stacked axes: each band its own ticks, in its own units
   const lb=R.left_band,rb=R.right_band,lo=Math.max(L[2],lb[0]),hi=Math.min(L[3],lb[1]);if(hi>lo)for(const t of nice(lo,hi))yTick(g,t,fmt(t),false);
   const a=Math.max(L[2],rb[0]),b=Math.min(L[3],rb[1]);if(b>a)for(const v of nice((a-R.offset)/R.scale,(b-R.offset)/R.scale))yTick(g,R.offset+R.scale*v,fmt(v),true);
   yTitle(g,f.y_label,lb[0],lb[1],false);yTitle(g,R.label,rb[0],rb[1],true);}
  else{for(const t of nice(L[2],L[3])){yTick(g,t,fmt(t),false);if(R){g.textAlign='left';g.fillText(fmt((t-R.offset)/R.scale),W()-m.r+5,sy(t)+4);}}
   yTitle(g,f.y_label,L[2],L[3],false);if(R)yTitle(g,R.label,L[2],L[3],true);}
  g.textAlign='center';g.fillText(f.x_label,(m.l+W()-m.r)/2,H()-8);
  g.save();g.beginPath();g.rect(m.l,m.t,W()-m.l-m.r,H()-m.t-m.b);g.clip();
  f.traces.forEach((t,k)=>{const col=t.color||P[k%P.length];if(t.heat){let lo=1/0,hi=-1/0;for(const v of t.heat)if(v!=null){lo=Math.min(lo,v);hi=Math.max(hi,v);}
    const cw=(W()-m.l-m.r)/t.cols,ch=(H()-m.t-m.b)/t.rows;for(let i=0;i<t.rows;i++)for(let j=0;j<t.cols;j++){g.fillStyle=viridis((t.heat[i*t.cols+j]-lo)/(hi-lo||1));g.fillRect(m.l+j*cw,H()-m.b-(i+1)*ch,cw+.5,ch+.5);}return;}
   g.strokeStyle=g.fillStyle=col;g.lineWidth=1.4;dash(g,t.style);
   if(t.marker){for(let i=0;i<t.x.length;i++)if(t.x[i]!=null&&t.y[i]!=null){g.beginPath();g.arc(sx(t.x[i]),sy(t.y[i]),2.5,0,7);g.fill();}}
   else{g.beginPath();let pen=false;for(let i=0;i<t.x.length;i++){if(t.x[i]==null||t.y[i]==null){pen=false;continue;}pen?g.lineTo(sx(t.x[i]),sy(t.y[i])):g.moveTo(sx(t.x[i]),sy(t.y[i]));pen=true;}g.stroke();}});
  g.restore();g.strokeStyle='#555';g.setLineDash([]);g.strokeRect(m.l,m.t,W()-m.l-m.r,H()-m.t-m.b);legend(g);}
 c.addEventListener('wheel',e=>{e.preventDefault();const k=e.deltaY>0?1.2:1/1.2,x=ix(e.offsetX),y=iy(e.offsetY);L=[x+(L[0]-x)*k,x+(L[1]-x)*k,y+(L[2]-y)*k,y+(L[3]-y)*k];draw();});
 let drag=null;c.addEventListener('mousedown',e=>drag=[e.offsetX,e.offsetY,L.slice()]);window.addEventListener('mouseup',()=>drag=null);
 c.addEventListener('mousemove',e=>{if(drag){const d=drag[2],dx=(e.offsetX-drag[0])/(W()-m.l-m.r)*(d[1]-d[0]),dy=(e.offsetY-drag[1])/(H()-m.t-m.b)*(d[3]-d[2]);L=[d[0]-dx,d[1]-dx,d[2]+dy,d[3]+dy];draw();return;}
  let best=null,bd=144;f.traces.forEach(t=>{if(!t.x)return;for(let i=0;i<t.x.length;i++){if(t.x[i]==null||t.y[i]==null)continue;const q=(sx(t.x[i])-e.offsetX)**2+(sy(t.y[i])-e.offsetY)**2;
   if(q<bd){bd=q;best=[t.label,t.x[i],t.right&&R?(t.y[i]-R.offset)/R.scale:t.y[i]];}}});
  o.textContent=best?(best[0]?best[0]+': ':'')+'('+fmt(best[1])+', '+fmt(best[2])+')':'x = '+fmt(ix(e.offsetX))+', y = '+fmt(iy(e.offsetY));});
 c.addEventListener('dblclick',()=>{L=home.slice();draw();});new ResizeObserver(draw).observe(c);return box;}
function show(p){root.replaceChildren();root.style.gridTemplateColumns='repeat('+Math.min(D.columns,Math.max(1,D.pages[p].figures.length))+',1fr)';for(const f of D.pages[p].figures)root.append(panel(f));[...tabs.children].forEach((b,i)=>b.classList.toggle('on',i===p));}
D.pages.forEach((p,i)=>{const b=document.createElement('button');b.textContent=p.label;b.onclick=()=>show(i);tabs.append(b);});if(D.pages.length<2)tabs.style.display='none';show(0);
)JS";

}  // namespace

std::string payloadHtml(const Payload& payload, const std::string& title) {
    auto data = json::dump(payloadJson(payload));
    std::string safe;  // keep "</script>" out of the embedded JSON
    for (const char c : data) { if (c == '<') safe += "\\u003c"; else safe += c; }
    std::string heading;
    for (const char c : title) { if (c == '<') heading += "&lt;"; else if (c == '&') heading += "&amp;"; else heading += c; }
    std::ostringstream out;
    out << "<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" "
           "content=\"width=device-width,initial-scale=1\"><title>" << heading << "</title><style>"
           "body{margin:0;font:14px Segoe UI,Helvetica,sans-serif;background:#f6f7f9;color:#16202b}"
           "h1{font-size:18px;margin:14px 18px 6px}#dsl-tabs{margin:0 18px}#dsl-tabs button{margin:0 4px 6px 0;"
           "padding:4px 10px;border:1px solid #c9d1da;background:#fff;border-radius:6px;cursor:pointer}"
           "#dsl-tabs button.on{background:#0f766e;color:#fff}#dsl-root{display:grid;gap:12px;padding:0 18px 18px}"
           ".fig{background:#fff;border:1px solid #c9d1da;border-radius:8px;position:relative}"
           ".fig canvas{width:100%;height:420px;display:block;cursor:crosshair}"
           ".readout{position:absolute;right:10px;bottom:6px;font-size:11px;color:#5b6876}</style></head><body>"
        << "<h1>" << heading << "</h1><div id=\"dsl-tabs\"></div><div id=\"dsl-root\"></div>"
        << "<script type=\"application/json\" id=\"dsl-data\">" << safe << "</script><script>" << kScript
        << "</script></body></html>";
    return out.str();
}

}  // namespace dsl
