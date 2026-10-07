
window.MathJax={tex:{inlineMath:[['$','$'],['\\(','\\)']],displayMath:[['$$','$$'],['\\[','\\]']]},options:{skipHtmlTags:['script','noscript','style','textarea','pre','code']}};


function fft(re,im){const n=re.length;
  for(let i=1,j=0;i<n;i++){let bit=n>>1;for(;j&bit;bit>>=1)j^=bit;j^=bit;
    if(i<j){let t=re[i];re[i]=re[j];re[j]=t;t=im[i];im[i]=im[j];im[j]=t;}}
  for(let len=2;len<=n;len<<=1){const ang=-2*Math.PI/len,half=len>>1;
    for(let i=0;i<n;i+=len)for(let k=0;k<half;k++){
      const wr=Math.cos(ang*k),wi=Math.sin(ang*k),ur=re[i+k],ui=im[i+k],
            ar=re[i+k+half],ai=im[i+k+half],vr=ar*wr-ai*wi,vi=ar*wi+ai*wr;
      re[i+k]=ur+vr;im[i+k]=ui+vi;re[i+k+half]=ur-vr;im[i+k+half]=ui-vi;}}}
function ifft(re,im){const n=re.length;for(let i=0;i<n;i++)im[i]=-im[i];fft(re,im);
  for(let i=0;i<n;i++){re[i]/=n;im[i]=-im[i]/n;}}
function analytic(x){const n=x.length,re=Float64Array.from(x),im=new Float64Array(n);
  fft(re,im);const h=new Float64Array(n);h[0]=1;
  if(n%2===0){h[n/2]=1;for(let i=1;i<n/2;i++)h[i]=2;}else{for(let i=1;i<=(n-1)/2;i++)h[i]=2;}
  for(let i=0;i<n;i++){re[i]*=h[i];im[i]*=h[i];}ifft(re,im);return{re:re,im:im};}
function mag(x){const n=x.length,re=Float64Array.from(x),im=new Float64Array(n);
  fft(re,im);const m=new Float64Array(n/2);
  for(let k=0;k<n/2;k++)m[k]=Math.hypot(re[k],im[k])*2/n;return m;}
function setupHiDPI(cv){
  const dpr=window.devicePixelRatio||1;
  if(!cv.dataset.baseH){cv.dataset.baseH=cv.getAttribute('height')||cv.clientHeight||300;}
  const h=parseInt(cv.dataset.baseH,10);
  const w=cv.clientWidth||cv.parentElement.clientWidth;
  cv.style.height=h+'px';
  const bw=Math.round(w*dpr),bh=Math.round(h*dpr);
  if(cv.width!==bw||cv.height!==bh){cv.width=bw;cv.height=bh;}
  const g=cv.getContext('2d');g.setTransform(dpr,0,0,dpr,0,0);
  return {g:g,W:w,H:h};
}
function poly(g,pts,color,width,dash){if(!pts.length)return;g.save();g.beginPath();
  if(dash)g.setLineDash(dash);
  for(let i=0;i<pts.length;i++){i?g.lineTo(pts[i][0],pts[i][1]):g.moveTo(pts[i][0],pts[i][1]);}
  g.strokeStyle=color;g.lineWidth=width||1.6;g.lineJoin='round';g.stroke();g.restore();}
function axis(g,x0,y0,x1,y1){poly(g,[[x0,y0],[x1,y1]],'#c9c4b6',1);}

/* ===== 图 8.1 SSB 频谱 ===== */
(function(){
  const cv=document.getElementById('c1'),sel=document.getElementById('sel1'),
        s=document.getElementById('s1'),v=document.getElementById('v1');
  const N=2048;
  function draw(){
    const {g,W,H}=setupHiDPI(cv);
    const fc=+s.value, mode=sel.value; v.textContent=fc;
    // 基带：不对称谱形（三根不等高谱线 + 形状）
    const m=new Float64Array(N);
    const comps=[[2,1.0],[4,0.72],[6,0.45],[8,0.22]];
    for(let i=0;i<N;i++){const t=i/N;let s0=0;
      for(const c of comps)s0+=c[1]*Math.cos(2*Math.PI*c[0]*t+c[0]*0.5);
      m[i]=s0;}
    const zm=analytic(m);
    const y=new Float64Array(N);
    for(let i=0;i<N;i++){const t=i/N;
      const c=Math.cos(2*Math.PI*fc*t), sn=Math.sin(2*Math.PI*fc*t);
      if(mode==='am') y[i]=m[i]*c;
      else if(mode==='usb') y[i]=m[i]*c-zm.im[i]*sn;
      else y[i]=m[i]*c+zm.im[i]*sn;}
    const M1=mag(m), M2=mag(y);
    let mx1=1e-9,mx2=1e-9;
    for(let k=0;k<50;k++)mx1=Math.max(mx1,M1[k]);
    for(let k=0;k<64;k++)mx2=Math.max(mx2,M2[k]);
    g.clearRect(0,0,W,H);
    const padL=14,padR=14,w=W-padL-padR;
    const KMAX=52, X=k=>padL+w*k/KMAX;
    const r1=H*0.33, a1=H*0.20, r2=H*0.90, a2=H*0.26;
    g.font='12px ui-sans-serif,system-ui,sans-serif';
    // 上：基带谱
    axis(g,padL,r1,W-padR,r1);
    for(let k=0;k<KMAX;k++){ if(M1[k]<mx1*0.012)continue;
      g.strokeStyle='#6b6f79';g.lineWidth=2.4;
      g.beginPath();g.moveTo(X(k),r1);g.lineTo(X(k),r1-a1*M1[k]/mx1);g.stroke();}
    g.fillStyle='#5d6472';g.fillText('基带 m(t) 的频谱（形状递减，便于分辨上下边带）',padL+2,r1-a1-8);
    g.fillText('0',padL-3,r1+16);
    // 下：调制后
    axis(g,padL,r2,W-padR,r2);
    const cl = mode==='am'?'#6b3fa8':(mode==='usb'?'#1f6f8b':'#b4531f');
    for(let k=0;k<KMAX;k++){ if(M2[k]<mx2*0.012)continue;
      g.strokeStyle=cl;g.lineWidth=2.4;
      g.beginPath();g.moveTo(X(k),r2);g.lineTo(X(k),r2-a2*M2[k]/mx2);g.stroke();}
    // 载波位置标记
    g.strokeStyle='#c0392b';g.setLineDash([4,3]);g.lineWidth=1.2;
    g.beginPath();g.moveTo(X(fc),r2+6);g.lineTo(X(fc),r2-a2-6);g.stroke();g.setLineDash([]);
    g.fillStyle='#c0392b';g.fillText('f_c',X(fc)-6,r2+17);
    g.fillStyle=cl;
    const nm = mode==='am'?'AM：上下边带都在，带宽 = 2B':
               (mode==='usb'?'SSB-USB：只保留 f_c 以上，带宽 = B':'SSB-LSB：只保留 f_c 以下，带宽 = B');
    g.fillText(nm,padL+2,r2-a2-8);
    g.fillStyle='#5d6472';g.fillText('频率 →',W-padR-44,r2+17);
  }
  sel.addEventListener('change',draw);s.addEventListener('input',draw);
  window.addEventListener('resize',draw);draw();
})();

/* ===== 图 8.2 轴承包络谱 ===== */
(function(){
  const cv=document.getElementById('c2');
  const s2=document.getElementById('s2'),v2=document.getElementById('v2');
  const s3=document.getElementById('s3'),v3=document.getElementById('v3');
  const s4=document.getElementById('s4'),v4=document.getElementById('v4');
  const N=4096, FS=20000;
  let seed=987654321;
  function rnd(){seed=(seed*1103515245+12345)&0x7fffffff;
    let u=seed/0x7fffffff;seed=(seed*1103515245+12345)&0x7fffffff;
    let u2=seed/0x7fffffff;
    return Math.sqrt(-2*Math.log(u+1e-12))*Math.cos(2*Math.PI*u2);}
  /* 简易 Butterworth 带通（两级双二阶级联，前后向滤波近似零相位） */
  function bandpass(x,lo,hi){
    const n=x.length;
    // 用 FFT 域理想带通（教学用途足够，且零相位）
    const re=Float64Array.from(x), im=new Float64Array(n);
    fft(re,im);
    for(let k=0;k<n;k++){
      const f=(k<=n/2? k : k-n)*FS/n;
      const af=Math.abs(f);
      let w=0;
      if(af>=lo&&af<=hi)w=1;
      else if(af>lo*0.7&&af<lo)w=0.5-0.5*Math.cos(Math.PI*(af-lo*0.7)/(lo*0.3));
      else if(af>hi&&af<hi*1.2)w=0.5+0.5*Math.cos(Math.PI*(af-hi)/(hi*0.2));
      re[k]*=w;im[k]*=w;
    }
    ifft(re,im); return re;
  }
  function draw(){
    const {g,W,H}=setupHiDPI(cv);
    const sev=(+s2.value)/100, nz=(+s3.value)/100, bpfo=+s4.value;
    v2.textContent=sev.toFixed(2);v3.textContent=nz.toFixed(2);v4.textContent=bpfo;
    seed=987654321;
    const x=new Float64Array(N);
    const fr=29.95, fres=3800, dec=900;   // 转频、共振频率、衰减率
    const T=1/bpfo;
    for(let i=0;i<N;i++){
      const t=i/FS;
      // 转频及其谐波（大幅低频干扰）
      let s0=1.0*Math.sin(2*Math.PI*fr*t)+0.45*Math.sin(2*Math.PI*2*fr*t+1.1)
             +0.30*Math.sin(2*Math.PI*3*fr*t+0.3)+0.55*Math.sin(2*Math.PI*9*fr*t+2.0);
      // 故障冲击：周期 T，激起衰减共振
      const k0=Math.floor(t/T);
      for(let k=Math.max(0,k0-2);k<=k0+1;k++){
        const dt=t-(k*T+0.01*T*Math.sin(k)); // 轻微随机滑移（抖动为周期的 1%）
        if(dt>=0&&dt<0.004)
          s0+=sev*1.15*Math.exp(-dec*dt)*Math.sin(2*Math.PI*fres*dt);
      }
      s0+=nz*0.55*rnd();
      x[i]=s0;
    }
    let mean=0;for(let i=0;i<N;i++)mean+=x[i];mean/=N;
    for(let i=0;i<N;i++)x[i]-=mean;
    // 带通 + 包络
    const xb=bandpass(x,2000,6000);
    const zb=analytic(xb);
    const env=new Float64Array(N);
    for(let i=0;i<N;i++)env[i]=Math.hypot(zb.re[i],zb.im[i]);
    let em=0;for(let i=0;i<N;i++)em+=env[i];em/=N;
    const env0=new Float64Array(N);
    for(let i=0;i<N;i++)env0[i]=(env[i]-em)*(0.5-0.5*Math.cos(2*Math.PI*i/(N-1)));
    const SP=mag(x), SE=mag(env0);
    g.clearRect(0,0,W,H);
    const padL=14,padR=14,w=W-padL-padR;
    g.font='11.5px ui-sans-serif,system-ui,sans-serif';
    const rows=[H*0.135,H*0.375,H*0.625,H*0.925];
    const amps=[H*0.075,H*0.085,H*0.075,H*0.105];
    // ① 原始波形
    axis(g,padL,rows[0],W-padR,rows[0]);
    let mx=1e-9;for(let i=0;i<N;i++)mx=Math.max(mx,Math.abs(x[i]));
    const p1=[];for(let i=0;i<N;i++)p1.push([padL+w*i/(N-1),rows[0]-amps[0]*x[i]/mx]);
    poly(g,p1,'#6b6f79',0.9);
    g.fillStyle='#5d6472';g.fillText('① 原始振动 x(t)（冲击被淹没）',padL+2,rows[0]-amps[0]-7);
    // ② 原始频谱（0~8000 Hz）
    const KM=Math.round(8000/FS*N);
    axis(g,padL,rows[1],W-padR,rows[1]);
    let sm=1e-9;for(let k=1;k<KM;k++)sm=Math.max(sm,SP[k]);
    // 共振带阴影
    const k1=Math.round(2000/FS*N),k2=Math.round(6000/FS*N);
    g.fillStyle='rgba(120,125,135,.13)';
    g.fillRect(padL+w*k1/KM,rows[1]-amps[1],w*(k2-k1)/KM,amps[1]);
    for(let k=1;k<KM;k++){const X=padL+w*k/KM;
      g.strokeStyle='#8a8f9a';g.lineWidth=0.8;
      g.beginPath();g.moveTo(X,rows[1]);g.lineTo(X,rows[1]-amps[1]*SP[k]/sm);g.stroke();}
    // BPFO 标记
    for(let h=1;h<=3;h++){const kk=h*bpfo/FS*N;
      if(kk<KM){g.strokeStyle='rgba(192,57,43,.55)';g.setLineDash([3,3]);g.lineWidth=1;
      g.beginPath();g.moveTo(padL+w*kk/KM,rows[1]);g.lineTo(padL+w*kk/KM,rows[1]-amps[1]);g.stroke();g.setLineDash([]);}}
    g.fillStyle='#5d6472';g.fillText('② x(t) 的频谱（0–8 kHz）　灰带 = 结构共振区　红虚线 = BPFO 谐波位置（看不见峰）',padL+2,rows[1]-amps[1]-7);
    // ③ 带通后 + 包络
    axis(g,padL,rows[2],W-padR,rows[2]);
    let bm=1e-9;for(let i=0;i<N;i++)bm=Math.max(bm,Math.abs(xb[i]));
    const p3=[],p3e=[];
    for(let i=0;i<N;i++){p3.push([padL+w*i/(N-1),rows[2]-amps[2]*xb[i]/bm]);
      p3e.push([padL+w*i/(N-1),rows[2]-amps[2]*env[i]/bm]);}
    poly(g,p3,'#9aa2ad',0.9); poly(g,p3e,'#c0392b',1.5);
    g.fillStyle='#c0392b';g.fillText('③ 带通 2–6 kHz 后的波形与包络（周期性冲击浮现）',padL+2,rows[2]-amps[2]-7);
    // ④ 包络谱（0~600 Hz）
    const KE=Math.round(600/FS*N);
    axis(g,padL,rows[3],W-padR,rows[3]);
    let em2=1e-9;for(let k=2;k<KE;k++)em2=Math.max(em2,SE[k]);
    for(let k=2;k<KE;k++){const X=padL+w*k/KE;
      g.strokeStyle='#1f6f8b';g.lineWidth=1.1;
      g.beginPath();g.moveTo(X,rows[3]);g.lineTo(X,rows[3]-amps[3]*SE[k]/em2);g.stroke();}
    for(let h=1;h<=4;h++){const f=h*bpfo;
      if(f<600){const X=padL+w*(f/FS*N)/KE;
      g.strokeStyle='rgba(192,57,43,.7)';g.setLineDash([4,3]);g.lineWidth=1.2;
      g.beginPath();g.moveTo(X,rows[3]+4);g.lineTo(X,rows[3]-amps[3]);g.stroke();g.setLineDash([]);
      g.fillStyle='#c0392b';g.font='10.5px ui-sans-serif,system-ui,sans-serif';
      g.fillText(h+'×',X-6,rows[3]+15);}}
    g.font='11.5px ui-sans-serif,system-ui,sans-serif';
    g.fillStyle='#1f6f8b';
    g.fillText('④ 包络谱（0–600 Hz）　红虚线 = BPFO = '+bpfo+' Hz 及谐波',padL+2,rows[3]-amps[3]-7);
    g.fillStyle='#5d6472';g.font='10.5px ui-sans-serif,system-ui,sans-serif';
    g.fillText('0',padL-2,H-2);g.fillText('600 Hz',W-padR-34,H-2);
  }
  s2.addEventListener('input',draw);s3.addEventListener('input',draw);s4.addEventListener('input',draw);
  window.addEventListener('resize',draw);draw();
})();
