// Which rainfall product where (/which-product/). Data: data/*.json, written by
// scripts/which_product_data.py and parked on blob by scripts/site_blob.py.
// fields.json (the per-cell map layers, ~1.7 MB gzipped) is the only large file and only the
// map section needs it, so the rest of the page is drawn from the small files first.
const DATA = ["summary", "outliers", "gauges", "gauge_viz", "gauges_3month", "tahmo_scores", "tahmo_funnel"];
const load = n => fetch(`data/${n}.json`).then(r => {
  if (!r.ok) throw new Error(`data/${n}.json: HTTP ${r.status}`);
  return r.json();
});
const fail = (where, what, e) => {
  const p = document.createElement("div"); p.className = "callout"; p.setAttribute("role", "alert");
  p.textContent = `${what} (${e.message}).`; where.prepend(p); throw e;
};
const fields = load("fields");
Promise.all(DATA.map(load)).then(d => {
  const initMap = main(...d);
  fields.then(initMap).catch(e => fail(document.getElementById("map"), "The map could not load or draw its data", e))
    .finally(() => document.getElementById("map-status").remove());
}).catch(e => fail(document.querySelector("main"), "The page could not load or draw its data; charts and tables below are missing", e));

// S summary, F per-cell fields, O outliers/PCA, G AMMA-CATCH scores, V AMMA-CATCH visuals,
// G3 AMMA-CATCH 3-month scores, T TAHMO regional scores, TV TAHMO selection funnel
function main(S, O, G, V, G3, T, TV) {
const P = S.products, LAB = S.labels, NRT = new Set(S.nrt), CORE = new Set(S.core), SCEN = Object.keys(S.scenarios);
const SHORT = {gpcc_full:"GPCC Full", gpcc_monitoring:"GPCC Monitoring", cru:"CRU TS", cpc:"CPC", precl:"PREC/L",
  chirp_v2:"CHIRP v2 (no stations)", chirp_v3:"CHIRP v3 (no stations)", chirps_v2:"CHIRPS v2", chirps_v3:"CHIRPS v3", imerg_late:"IMERG Late", imerg_final:"IMERG Final", era5:"ERA5"};
const SCEN_SHORT = {sat: "A: satellites share mistakes", smooth: "B: GPCC and ERA5 share mistakes"};
const css = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
const dec = (s, T) => { const b = atob(s), u = new Uint8Array(b.length); for (let i = 0; i < b.length; i++) u[i] = b.charCodeAt(i); return new T(u.buffer); };
const farP = dec(O.cell.farthest, Int8Array), farD = dec(O.cell.far_dist, Uint8Array), disg = dec(O.cell.disagree, Uint8Array);
const GCLS = ["never gauged", "gauged <50% of months", "gauged ≥50% of months"];
const SCALES = ["monthly", "3-month"];
const fmt = (v, d = 2) => v == null || Number.isNaN(v) ? "–" : v.toFixed(d);
const pct = v => Math.round(v * 100) + "%";
let scen = S.primary;

const hex = h => { h = h.replace("#", ""); return [0, 2, 4].map(i => parseInt(h.slice(i, i + 2), 16)); };
const mix = (a, b, t) => a.map((x, i) => Math.round(x + (b[i] - x) * t));
const rgb = c => `rgb(${c[0]},${c[1]},${c[2]})`;
const SEQ = ["#f2f7d9", "#a8d8a0", "#3fae8f", "#1c6f8a", "#2a2f6b"].map(hex);
function seq(t) { t = Math.max(0, Math.min(1, t)); const k = t * (SEQ.length - 1), i = Math.min(Math.floor(k), SEQ.length - 2); return mix(SEQ[i], SEQ[i + 1], k - i); }
function div(t) { const n = hex(css("--neg")), m = hex(css("--mid")), p = hex(css("--pos")); t = Math.max(-1, Math.min(1, t)); return t < 0 ? mix(m, n, -t) : mix(m, p, t); }
const REG_COL = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#1f78b4", "#b15928", "#6a3d9a"];
const PROD_COL = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02"];
const GAUGE_COL = ["#e6ab02", "#7fb8a4", "#1e795f"];
const inkFor = c => (0.299 * c[0] + 0.587 * c[1] + 0.114 * c[2]) > 150 ? "#1f2324" : "#ffffff";

function seg(el, opts, cur, on, labels) {
  el.innerHTML = ""; for (const o of opts) { const b = document.createElement("button"); b.type = "button"; b.textContent = labels ? labels[o] : o;
    b.setAttribute("aria-pressed", o === cur); b.onclick = () => { [...el.children].forEach(c => c.setAttribute("aria-pressed", c === b)); on(o); }; el.appendChild(b); }
}
function cell(v, color, d = 2) { const td = document.createElement("td"); if (v == null || Number.isNaN(v)) { td.textContent = "–"; td.style.color = "var(--n7)"; return td; }
  td.textContent = fmt(v, d); td.style.background = rgb(color); td.style.color = inkFor(color); return td; }
const th = (txt, vert, h = 130) => vert ? `<th style="writing-mode:vertical-rl;transform:rotate(180deg);height:${h}px">${txt}</th>` : `<th>${txt}</th>`;

// ---- static numbers in prose
document.getElementById("fill-n").textContent = S.chirps_v3_unscored_cells.toLocaleString();
document.getElementById("sil").textContent = S.k_silhouette[S.k].toFixed(2);
document.getElementById("px-stable").textContent = pct(S.pixel_best_nrt_stable_frac);
(function () {
  const t = document.getElementById("tetrad-table");
  t.innerHTML = "<tr><th>Cells</th>" + GCLS.map(g => `<th>${g}</th>`).join("") + "<th>Trend</th></tr>";
  const row = (lab, vals) => { const tr = document.createElement("tr"); tr.innerHTML = `<td class="rowh">${lab}</td>` + vals.map(v => `<td class="num">${v == null ? "–" : (v >= 0 ? "+" : "") + v.toFixed(3)}</td>`).join("") +
    `<td style="text-align:left">${vals[2] > vals[0] ? "rises → satellites share mistakes" : "falls → GPCC and ERA5 do"}</td>`; t.appendChild(tr); };
  row("Inconsistency, all land", GCLS.map(g => S.tetrad[g].delta));
  row("Inconsistency, Africa", GCLS.map(g => S.tetrad_africa[g].delta));
  const c = document.createElement("tr"); c.innerHTML = `<td class="rowh">Control (should be ≈ 0)</td>` + GCLS.map(g => `<td class="num">${S.tetrad[g].control.toFixed(3)}</td>`).join("") + "<td></td>"; t.appendChild(c);
})();

// ---- inputs table + diagram (scenario-aware "scored against")
function drawInputs() {
  const t = document.getElementById("inputs-table"), fams = ["G", "C", "I", "E"], used = S.scenarios[scen].used;
  t.innerHTML = "<tr><th>Product</th>" + fams.map(f => `<th>${f}</th>`).join("") + "<th>Latency</th><th>Scored against</th></tr>";
  for (const p of P) { const tr = document.createElement("tr"); const h = document.createElement("td"); h.className = "rowh" + (CORE.has(p) ? " core" : ""); h.textContent = SHORT[p]; tr.appendChild(h);
    for (const f of fams) { const td = document.createElement("td"); td.textContent = S.inputs[p].includes(f) ? "■" : ""; td.style.textAlign = "center"; td.style.color = "var(--b6)"; tr.appendChild(td); }
    const lat = document.createElement("td"); lat.textContent = S.latency[p]; lat.style.textAlign = "left"; tr.appendChild(lat);
    const u = used[p]; const part = CORE.has(p) ? u.map(tr_ => tr_.map(k => SHORT[k]).join(" + ")).join("; ") : u.map(k => SHORT[k]).join(", ");
    const pt = document.createElement("td"); pt.textContent = part; pt.style.textAlign = "left"; tr.appendChild(pt); t.appendChild(tr); }
}
(function () {
  const g = document.getElementById("core-boxes"); const xs = [10, 102, 194, 286], cx = i => xs[i] + 42;
  S.core.forEach((k, i) => { const x = xs[i];
    g.insertAdjacentHTML("beforeend", `<line x1="190" y1="52" x2="${cx(i)}" y2="150" stroke="var(--b5)" stroke-width="1.5"/>
      <rect x="${x}" y="150" width="84" height="38" rx="4" fill="var(--surface)" stroke="var(--n3)"/>
      <text x="${cx(i)}" y="173" text-anchor="middle" font-size="11" font-weight="500">${SHORT[k]}</text>`); });
  // candidate shared errors: CHIRP v3 (no stations)–IMERG Late (idx 1,2) and GPCC–ERA5 (idx 0,3)
  g.insertAdjacentHTML("beforeend", `<path d="M ${cx(1)} 188 Q ${(cx(1) + cx(2)) / 2} 222 ${cx(2)} 188" fill="none" stroke="var(--neg)" stroke-width="1.5" stroke-dasharray="4 3"/>
    <text x="${(cx(1) + cx(2)) / 2}" y="222" text-anchor="middle" font-size="10">shared infrared mistakes</text>`);
})();

// ---- dependency matrices (residuals: primary scenario)
const depRegions = Object.keys(S.dependency.monthly);
let depScale = "monthly";
const depSel = document.getElementById("dep-region");
depRegions.forEach(r => depSel.add(new Option(`${r} (${S.dependency.monthly[r].cells.toLocaleString()} cells)`, r)));
depSel.value = "never gauged";
function matrix(tbl, M, order, colorFn, cap) {
  tbl.innerHTML = `<caption class="muted" style="text-align:left;caption-side:top">${cap}</caption>`;
  const hr = document.createElement("tr"); hr.innerHTML = "<th></th>" + order.map(p => th(SHORT[p], true, 92)).join(""); tbl.appendChild(hr);
  for (const a of order) { const tr = document.createElement("tr"); const h = document.createElement("td"); h.className = "rowh"; h.textContent = SHORT[a]; tr.appendChild(h);
    for (const b of order) { const v = M[P.indexOf(a)][P.indexOf(b)]; tr.appendChild(v == null ? document.createElement("td") : cell(v, colorFn(v))); } tbl.appendChild(tr); }
}
function drawDep() {
  const D = S.dependency[depScale][depSel.value];
  matrix(document.getElementById("dep-corr"), D.corr, D.order, v => seq(v), "Anomaly correlation, ordered by hierarchical clustering");
  matrix(document.getElementById("dep-resid"), D.resid, D.order, v => div(v / 0.4), "Residual after the estimated common signal (blue = shared error)");
}
seg(document.getElementById("dep-scale"), SCALES, depScale, s => { depScale = s; drawDep(); });
depSel.onchange = drawDep;

// ---- skill table
let skillScale = "monthly";
function drawSkill() {
  const T = S.skill[scen][skillScale], regs = Object.keys(T), t = document.getElementById("skill-table");
  t.innerHTML = "<tr><th>Product</th>" + regs.map(r => th(`${r} <span class="muted">(${S.regions[r].toLocaleString()})</span>`, true, 160)).join("") + "</tr>";
  const order = [...P].sort((a, b) => (T.All[b].median ?? -1) - (T.All[a].median ?? -1));
  for (const p of order) { const tr = document.createElement("tr"); const h = document.createElement("td"); h.className = "rowh" + (CORE.has(p) ? " core" : "");
    h.textContent = SHORT[p]; if (NRT.has(p)) h.style.color = "var(--pos)"; tr.appendChild(h);
    for (const r of regs) { const v = T[r][p].median; tr.appendChild(cell(v, seq(v ?? 0))); } t.appendChild(tr); }
}
seg(document.getElementById("skill-scale"), SCALES, skillScale, s => { skillScale = s; drawSkill(); });

// ---- decision table: best NRT product per region under each scenario
let decScale = "monthly";
function drawDecision() {
  const T = S.skill[scen][decScale], regs = Object.keys(T), t = document.getElementById("dec-table");
  t.innerHTML = "<tr><th>Region</th><th>Best near-real-time</th><th>Runner-up</th><th>Lead</th><th>Best gauge analysis available within ~2 months</th></tr>";
  for (const r of regs) {
    const o = S.nrt.filter(p => T[r][p].median != null).sort((a, b) => T[r][b].median - T[r][a].median);
    const lead = T[r][o[0]].median - T[r][o[1]].median;
    const gBest = ["gpcc_monitoring", "cpc"].sort((a, b) => T[r][b].median - T[r][a].median)[0];
    const tr = document.createElement("tr");
    tr.innerHTML = `<td class="rowh">${r}</td><td style="text-align:left"><b>${SHORT[o[0]]}</b> <span class="num muted">${fmt(T[r][o[0]].median)}</span></td>` +
      `<td style="text-align:left">${SHORT[o[1]]} <span class="num muted">${fmt(T[r][o[1]].median)}</span></td>` +
      `<td style="text-align:left"><span class="chip ${lead < 0.02 ? "warn" : "ok"}">${lead < 0.02 ? "tie" : "+" + fmt(lead)}</span></td>` +
      `<td style="text-align:left">${SHORT[gBest]} <span class="num muted">${fmt(T[r][gBest].median)}</span></td>`;
    t.appendChild(tr);
  }
  document.getElementById("dec-note").textContent = "Median ρ² over the region's cells. Gauge analyses are shown for comparison: GPCC Monitoring arrives about two months late, CPC within a day.";
}
seg(document.getElementById("dec-scale"), SCALES, decScale, s => { decScale = s; drawDecision(); });

// ---- map (drawn once fields.json arrives; drawMap is a no-op until then)
let drawMap = () => {};
function initMap(F) {
const iy = dec(F.iy, Int16Array), ix = dec(F.ix, Int16Array), regime = dec(F.regime, Int8Array), bestNrt = dec(F.best_nrt, Int8Array),
  bestAlt = dec(F.best_alt, Int8Array), robust = dec(F.robust, Int8Array),
  margin = dec(F.margin, Uint8Array), pc1 = dec(F.pc1, Uint8Array), pc1c = dec(F.pc1_core, Uint8Array), gauge = dec(F.gauge, Int8Array);
const rho = {}; for (const sc in F.rho2) { rho[sc] = {}; for (const s in F.rho2[sc]) { rho[sc][s] = {}; for (const p in F.rho2[sc][s]) rho[sc][s][p] = dec(F.rho2[sc][s][p], Uint8Array); } }
const cv = document.getElementById("map-canvas"), ctx = cv.getContext("2d");
let y0 = Infinity, y1 = -Infinity; for (let i = 0; i < F.n; i++) { y0 = Math.min(y0, iy[i]); y1 = Math.max(y1, iy[i]); }
const H = y1 - y0 + 1; cv.height = H;
const grid = new Int32Array(720 * H).fill(-1); for (let i = 0; i < F.n; i++) grid[(y1 - iy[i]) * 720 + ix[i]] = i;
const LAYERS = [["best", "Best near-real-time product"], ["margin", "Lead of best over second (ρ²)"],
  ["far", "Product farthest from the other eleven"], ["fardist", "How far the farthest product is (1 − r)"], ["disagree", "Average disagreement between products"],
  ["pc1", "Consensus, all 12 products (PC1 share)"], ["pc1c", "Consensus, 4 core products (PC1 share)"], ["gauge", "GPCC gauge presence"], ["regime", "Exploratory regimes (k-means)"],
  ...P.map(p => ["rho:" + p, "ρ² " + SHORT[p] + ""])];
const ls = document.getElementById("map-layer"); LAYERS.forEach(([k, l]) => ls.add(new Option(l, k)));
let mapScale = "monthly";
const nrtList = S.nrt;
function layerColor(k, i) {
  if (k === "regime") return regime[i] < 0 ? null : hex(REG_COL[regime[i]]);
  if (k === "far") return farP[i] < 0 ? null : hex(PCOL[P[farP[i]]]);
  if (k === "fardist") return farD[i] === 255 ? null : seq(farD[i] / 254 / 0.6);
  if (k === "disagree") return disg[i] === 255 ? null : seq((disg[i] / 254 - 0.3) / 0.5);
  if (k === "best") return bestNrt[i] < 0 ? null : hex(PROD_COL[bestNrt[i]]);
  if (k === "bestalt") return bestAlt[i] < 0 ? null : hex(PROD_COL[bestAlt[i]]);
  if (k === "robust") return robust[i] ? hex(PROD_COL[bestNrt[i]]) : hex(css("--n2"));
  if (k === "gauge") return gauge[i] < 0 ? null : hex(GAUGE_COL[gauge[i]]);
  const u = k === "margin" ? margin[i] : k === "pc1" ? pc1[i] : k === "pc1c" ? pc1c[i] : rho[scen][mapScale][k.slice(4)][i];
  if (u === 255) return null; const v = u / 254;
  return seq(k.startsWith("pc1") ? (v - 0.25) / 0.75 : v);
}
function legend(k) {
  const el = document.getElementById("map-legend");
  const ramp = (lo, hi, lab) => `<span>${lab}</span><span class="num">${lo}</span><span class="ramp" style="background:linear-gradient(90deg,${[0, .25, .5, .75, 1].map(t => rgb(seq(t))).join(",")})"></span><span class="num">${hi}</span>`;
  const prods = nrtList.map((p, i) => `<span><span class="sw" style="background:${PROD_COL[i]}"></span>${SHORT[p]}</span>`).join("");
  if (k === "regime") el.innerHTML = S.regimes.map(r => `<span><span class="sw" style="background:${REG_COL[r.id]}"></span>Regime ${r.id + 1} (${pct(r.area_frac)})</span>`).join("");
  else if (k === "robust") el.innerHTML = prods + `<span><span class="sw" style="background:${css("--n2")};border:1px solid var(--n3)"></span>A and B disagree</span>`;
  else if (k === "best" || k === "bestalt") el.innerHTML = prods;
  else if (k === "far") el.innerHTML = P.map(p => `<span><span class="sw" style="background:${PCOL[p]}"></span>${SHORT[p]}</span>`).join("");
  else if (k === "fardist") el.innerHTML = ramp("0", "0.6+", "1 − correlation with the others' average");
  else if (k === "disagree") el.innerHTML = ramp("0.3", "0.8+", "SD across products (standardised)");
  else if (k === "gauge") el.innerHTML = GCLS.map((g, i) => `<span><span class="sw" style="background:${GAUGE_COL[i]}"></span>${g}</span>`).join("");
  else if (k === "margin") el.innerHTML = ramp("0", "0.5", "ρ² lead") + `<span class="muted">pale = near-tie</span>`;
  else if (k.startsWith("pc1")) el.innerHTML = ramp("0.25", "1", "share of variance on PC1");
  else el.innerHTML = ramp("0", "1", "ρ²");
  el.insertAdjacentHTML("beforeend", `<span><span class="sw" style="background:${css("--n3")}"></span>not scored</span>`);
}
drawMap = function () {
  const k = ls.value, img = ctx.createImageData(720, H), oc = hex(css("--ocean")), na = hex(css("--n3"));
  for (let p = 0; p < 720 * H; p++) { const i = grid[p]; const c = i < 0 ? oc : (layerColor(k, i) || na);
    img.data[p * 4] = c[0]; img.data[p * 4 + 1] = c[1]; img.data[p * 4 + 2] = c[2]; img.data[p * 4 + 3] = 255; }
  ctx.putImageData(img, 0, 0); legend(k);
};
ls.onchange = drawMap; seg(document.getElementById("map-scale"), SCALES, mapScale, s => { mapScale = s; drawMap(); });
let lastCell = -1;
function drawPanel(i) {
  const pn = document.getElementById("cell-panel"); lastCell = i;
  if (i == null || i < 0) { pn.innerHTML = "<h4>No data here</h4><p class='muted'>Ocean, or a cell outside the common domain.</p>"; return; }
  const lat = -89.75 + 0.5 * iy[i], lon = -179.75 + 0.5 * ix[i];
  const rows = [scen].map(sc => { const R = rho[sc][mapScale]; const ord = [...P].sort((a, b) => (R[b][i] === 255 ? -1 : R[b][i]) - (R[a][i] === 255 ? -1 : R[a][i]));
    return `<h4 style="margin-top:10px">${mapScale}</h4><div class="bars">` + ord.map(p => { const u = R[p][i], v = u === 255 ? null : u / 254;
      return `<span class="lab">${SHORT[p]}</span><span><span class="bar${NRT.has(p) ? " nrt" : ""}" style="display:block;width:${v == null ? 0 : v * 100}%"></span></span><span class="v">${fmt(v)}</span>`; }).join("") + "</div>"; }).join("");
  pn.innerHTML = `<h4>${lat.toFixed(2)}°, ${lon.toFixed(2)}°</h4><div class="kv"><dt>Gauges</dt><dd>${GCLS[gauge[i]] || "–"}</dd>
    <dt>Best for triggers</dt><dd>${SHORT[nrtList[bestNrt[i]]]}</dd><dt>Lead</dt><dd class="num">${fmt(margin[i] === 255 ? null : margin[i] / 254 * 0.5)}</dd></div>${rows}
    <p class="muted">Blue bars: near-real-time products. One cell is noisy; compare with its region.</p>`;
}
cv.addEventListener("click", e => { const r = cv.getBoundingClientRect(), x = Math.floor((e.clientX - r.left) / r.width * 720), y = Math.floor((e.clientY - r.top) / r.height * H);
  drawPanel(grid[y * 720 + x]); });
drawMap();
}

// ---- robustness
(function () {
  const sel = document.getElementById("rob-region"); ["All", "never gauged", "gauged ≥50% of months", "Africa"].forEach(r => sel.add(new Option(r, r)));
  const t = document.getElementById("rob-table");
  function draw() { const reg = sel.value, alts = Object.keys(S.sensitivity.monthly);
    t.innerHTML = "<tr><th>Product</th><th>Main result</th><th>What-if: GPCC and ERA5 share mistakes</th>" + alts.map(a => `<th>${a}</th>`).join("") + "</tr>";
    for (const p of P) { const tr = document.createElement("tr"); const h = document.createElement("td"); h.className = "rowh"; h.textContent = SHORT[p]; tr.appendChild(h);
      for (const sc of SCEN) { const v = S.skill[sc].monthly[reg][p].median; tr.appendChild(cell(v, seq(v ?? 0))); }
      for (const a of alts) { const v = S.sensitivity.monthly[a][reg][p]; tr.appendChild(cell(v, seq(v ?? 0))); } t.appendChild(tr); } }
  sel.onchange = draw; draw();
  const cl = P.filter(p => S.clipped_high.monthly[p] > 0.02);
  document.getElementById("rob-note").textContent = `Monthly, median ρ² over cells. ${cl.length ? "Estimates above 1 were clipped in: " + cl.map(p => `${SHORT[p]} ${pct(S.clipped_high.monthly[p])}`).join(", ") + " of cells." : "Under 2% of estimates exceeded 1 and were clipped for every product."}`;
})();

// ---- PCA biplot, outlier shares, disagreement series
const PCOL = {gpcc_full:"#1b7f5b", gpcc_monitoring:"#3fae8f", cru:"#7fc8a9", cpc:"#2f6f6a", precl:"#9bc4a0",
  chirp_v2:"#f08c2e", chirp_v3:"#d0661b", chirps_v2:"#b5432a", chirps_v3:"#8a2b1d", imerg_late:"#8e6cc4", imerg_final:"#5a3d99", era5:"#2a78d6"};
const pcaSel = document.getElementById("pca-region");
Object.keys(O.pca).forEach(r => pcaSel.add(new Option(`${r} (${O.pca[r].cells.toLocaleString()} cells)`, r)));
let pcaAxes = "PC1 × PC2", pcaMode = "Raw";
const MODE_NOTE = {
  "Raw": "Raw: every product is mostly the same rain signal, so PC1 takes most of the variance, all products crowd to the right, and they fall on an arc near the unit circle (each product's squared loadings over all components sum to 1). This is the same geometry PCA gives on closed (compositional) data without a log-ratio transform: one dominant shared component sets the shape.",
  "Deviations from consensus": "Deviations from consensus: each product's anomaly minus the average of all twelve in that cell and month, the analogue of a centred log-ratio for closed data. The shared rain signal is gone, so the plot shows only how products differ, and distance from the centre is how much a product departs from the crowd. Deviations sum to zero across products, which forces one component to zero and makes loadings contrasts between products."};
function drawBiplot() {
  const reg = pcaSel.value, pc = (pcaMode === "Raw" ? O.pca : O.pca_dev)[reg], [a, b] = pcaAxes === "PC1 × PC2" ? [0, 1] : [1, 2];
  const svg = document.getElementById("biplot"), W = 520, H = 440, m = 46;
  const xs = P.map(p => pc.loadings[p][a]), ys = P.map(p => pc.loadings[p][b]);
  const pad = 0.08, x0 = Math.min(0, ...xs) - pad, x1 = Math.max(0, ...xs) + pad, y0 = Math.min(0, ...ys) - pad, y1 = Math.max(0, ...ys) + pad;
  const sx = v => m + (v - x0) / (x1 - x0) * (W - 2 * m), sy = v => H - m - (v - y0) / (y1 - y0) * (H - 2 * m);
  const tick = (v0, v1) => { const out = []; for (let v = Math.ceil(v0 * 5) / 5; v <= v1 + 1e-9; v += 0.2) out.push(+v.toFixed(1)); return out; };
  let g = `<rect x="${m}" y="${m}" width="${W - 2 * m}" height="${H - 2 * m}" fill="var(--surface)" stroke="var(--n2)"/>`;
  for (const v of tick(x0, x1)) g += `<line x1="${sx(v)}" x2="${sx(v)}" y1="${m}" y2="${H - m}" stroke="var(--n1)"/><text x="${sx(v)}" y="${H - m + 14}" font-size="10" text-anchor="middle">${v}</text>`;
  for (const v of tick(y0, y1)) g += `<line x1="${m}" x2="${W - m}" y1="${sy(v)}" y2="${sy(v)}" stroke="var(--n1)"/><text x="${m - 6}" y="${sy(v) + 3}" font-size="10" text-anchor="end">${v}</text>`;
  g += `<line x1="${sx(0)}" x2="${sx(0)}" y1="${m}" y2="${H - m}" stroke="var(--n3)"/><line x1="${m}" x2="${W - m}" y1="${sy(0)}" y2="${sy(0)}" stroke="var(--n3)"/>`;
  g += `<text x="${W / 2}" y="${H - 8}" font-size="11" text-anchor="middle">PC${a + 1} (${Math.round(pc.var_frac[a] * 100)}% of variance)</text>`;
  g += `<text transform="translate(12 ${H / 2}) rotate(-90)" font-size="11" text-anchor="middle">PC${b + 1} (${Math.round(pc.var_frac[b] * 100)}% of variance)</text>`;
  // labels: nudge to reduce overlap
  const pts = P.map((p, i) => ({p, x: sx(xs[i]), y: sy(ys[i]), ly: sy(ys[i])})).sort((u, v) => u.ly - v.ly);
  for (let i = 1; i < pts.length; i++) for (let j = 0; j < i; j++) if (Math.abs(pts[i].x - pts[j].x) < 70 && Math.abs(pts[i].ly - pts[j].ly) < 12) pts[i].ly = pts[j].ly + 12;
  for (const q of pts) {
    const left = q.x > W - m - 80;
    g += `<line x1="${sx(0)}" y1="${sy(0)}" x2="${q.x}" y2="${q.y}" stroke="${PCOL[q.p]}" stroke-opacity=".35"/>`;
    g += `<circle cx="${q.x}" cy="${q.y}" r="5" fill="${PCOL[q.p]}" stroke="var(--surface)"/>`;
    g += `<text x="${q.x + (left ? -8 : 8)}" y="${q.ly + 4}" font-size="11" text-anchor="${left ? "end" : "start"}" font-weight="500">${SHORT[q.p]}</text>`;
  }
  svg.innerHTML = g;
  // outlier shares
  const sh = O.share_farthest[reg], md = O.median_distance[reg];
  const ord = [...P].sort((u, v) => sh[v] - sh[u]);
  document.getElementById("far-bars").innerHTML = ord.map(p => `<span class="lab">${SHORT[p]}</span><span><span class="bar" style="display:block;background:${PCOL[p]};width:${sh[p] * 100}%"></span></span><span class="v">${pct(sh[p])}</span>`).join("");
  document.getElementById("pca-note").textContent = `Distance from the others = 1 − correlation with the average of the other eleven. Median distance in this region: ` +
    [...P].sort((u, v) => md[v] - md[u]).slice(0, 4).map(p => `${SHORT[p]} ${fmt(md[p])}`).join(", ") + `.` + (pcaMode === "Raw" ? " PC1 is the shared rain signal; PC2 separates gauge analyses (positive) from satellite products without stations (negative)." : " PC1 of the deviations is the main axis of disagreement: gauge analyses against CHIRP, with CHIRPS between.");
  drawSpread(reg);
}
function drawSpread(reg) {
  const svg = document.getElementById("spread-ts"), W = 760, H = 180, m = {l: 40, r: 10, t: 10, b: 24};
  const ys = O.spread_ts.series[reg] || O.spread_ts.series.All, ts = O.spread_ts.times, n = ys.length;
  const vals = ys.filter(v => v != null), lo = Math.floor(Math.min(...vals) * 10) / 10, hi = Math.ceil(Math.max(...vals) * 10) / 10;
  const sx = i => m.l + i / (n - 1) * (W - m.l - m.r), sy = v => H - m.b - (v - lo) / (hi - lo) * (H - m.t - m.b);
  let g = "";
  for (let v = lo; v <= hi + 1e-9; v += 0.1) g += `<line x1="${m.l}" x2="${W - m.r}" y1="${sy(v)}" y2="${sy(v)}" stroke="var(--n1)"/><text x="${m.l - 6}" y="${sy(v) + 3}" font-size="10" text-anchor="end">${v.toFixed(1)}</text>`;
  ts.forEach((t, i) => { if (t.endsWith("-01") && +t.slice(0, 4) % 2 === 1) g += `<text x="${sx(i)}" y="${H - 6}" font-size="10" text-anchor="middle">${t.slice(0, 4)}</text>`; });
  let d = "", area = ""; ys.forEach((v, i) => { if (v == null) return; d += (d ? "L" : "M") + sx(i).toFixed(1) + " " + sy(v).toFixed(1); });
  area = d + `L${sx(n - 1)} ${H - m.b}L${sx(0)} ${H - m.b}Z`;
  g += `<path d="${area}" fill="var(--b05)"/><path d="${d}" fill="none" stroke="var(--b6)" stroke-width="1.3"/>`;
  const top = ys.map((v, i) => [v ?? -1, i]).sort((u, v) => v[0] - u[0]).slice(0, 3);
  for (const [v, i] of top) g += `<circle cx="${sx(i)}" cy="${sy(v)}" r="3" fill="var(--neg)"/><text x="${sx(i)}" y="${sy(v) - 6}" font-size="10" text-anchor="middle">${ts[i]}</text>`;
  svg.innerHTML = g;
}
pcaSel.onchange = drawBiplot;
seg(document.getElementById("pca-mode"), ["Raw", "Deviations from consensus"], pcaMode, s => { pcaMode = s; document.getElementById("pca-mode-note").textContent = MODE_NOTE[s]; drawBiplot(); });
document.getElementById("pca-mode-note").textContent = MODE_NOTE[pcaMode];
seg(document.getElementById("pca-axes"), ["PC1 × PC2", "PC2 × PC3"], pcaAxes, s => { pcaAxes = s; drawBiplot(); });
drawBiplot();

// ---- gauge validation tables
(function () {
  const t = document.getElementById("gauge-table"), ps = Object.keys(G.products).sort((a, b) => G.products[b].r2_vs_gauge - G.products[a].r2_vs_gauge);
  const ci = c => `<span class="muted num" style="font-size:11px"> ${fmt(c[0])}–${fmt(c[1])}</span>`;
  t.innerHTML = "<tr><th>Product</th><th>r²</th><th>Spearman</th><th>Dry-month hit rate</th><th>False alarms</th><th>Available</th></tr>";
  for (const p of ps) { const v = G.products[p], tr = document.createElement("tr");
    const h = document.createElement("td"); h.className = "rowh"; h.textContent = SHORT[p]; if (NRT.has(p)) h.style.color = "var(--pos)"; tr.appendChild(h);
    for (const [k, sc] of [["r2_vs_gauge", 0.7], ["spearman", 0.85], ["dry_hit", 0.8]]) { const td = cell(v[k], seq(v[k] / sc)); td.innerHTML += ci(v[k === "r2_vs_gauge" ? "r2_ci" : k + "_ci"]); tr.appendChild(td); }
    const fa = document.createElement("td"); fa.className = "num"; fa.textContent = pct(v.dry_false_alarm); tr.appendChild(fa);
    const la = document.createElement("td"); la.style.textAlign = "left"; la.textContent = S.latency[p]; tr.appendChild(la);
    t.appendChild(tr); }
  const m = document.getElementById("gauge-vs-model");
  m.innerHTML = "<tr><th>Product</th><th>Gauge-based ρ²</th><th>Estimate: if satellites share mistakes</th><th>Estimate: if GPCC and ERA5 share mistakes</th></tr>";
  for (const p of ps) { const v = G.products[p], tr = document.createElement("tr");
    const h = document.createElement("td"); h.className = "rowh"; h.textContent = SHORT[p]; tr.appendChild(h);
    tr.appendChild(cell(v.gauge_based, seq(v.gauge_based)));
    for (const k of ["model_sat", "model_smooth"]) { const td = cell(v[k], seq(v[k])); const d = v[k] - v.gauge_based;
      td.innerHTML += `<span class="num" style="font-size:11px;opacity:.8"> (${d >= 0 ? "+" : ""}${fmt(d)})</span>`; tr.appendChild(td); }
    m.appendChild(tr); }
})();

// ---- gauge results by month and period
(function () {
  const MON = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"];
  const order = Object.keys(G.products).sort((a, b) => G.products[b].r2_vs_gauge - G.products[a].r2_vs_gauge);
  const sites = Object.keys(G.by_month).filter(s => Object.values(G.by_month[s]).some(g => g));
  const gs = document.getElementById("gm-site"), ps = document.getElementById("gp-site");
  sites.forEach(s => { gs.add(new Option(s === "All" ? "All gauged cells" : s, s)); ps.add(new Option(s === "All" ? "All gauged cells" : s, s)); });
  function drawMonth() {
    const bm = G.by_month[gs.value], months = Object.keys(bm).filter(m => bm[m]).map(Number);
    const t = document.getElementById("gauge-month");
    t.innerHTML = "<tr><th>Product</th>" + months.map(m => `<th>${MON[m - 1]}<br><span class="muted" style="font-weight:400">n=${bm[m].n}</span></th>`).join("") + "</tr>";
    const best = {}; for (const m of months) best[m] = S.nrt.reduce((a, b) => bm[m].products[b].r2 > bm[m].products[a].r2 ? b : a);
    for (const p of order) { const tr = document.createElement("tr"); const h = document.createElement("td"); h.className = "rowh"; h.textContent = SHORT[p]; if (NRT.has(p)) h.style.color = "var(--pos)"; tr.appendChild(h);
      for (const m of months) { const v = bm[m].products[p], td = cell(v.r2, seq(v.r2 / 0.8)); td.title = `${SHORT[p]}, ${MON[m - 1]}: r² ${fmt(v.r2)} (90% range ${fmt(v.lo)}–${fmt(v.hi)})`;
        if (best[m] === p) td.style.outline = "2px solid var(--n9)"; tr.appendChild(td); }
      t.appendChild(tr); }
  }
  function drawPeriod() {
    const bp = G.by_period[ps.value], ws = Object.keys(bp).filter(w => bp[w]);
    const t = document.getElementById("gauge-period");
    t.innerHTML = "<tr><th>Product</th>" + ws.map(w => `<th>${w}<br><span class="muted" style="font-weight:400">n=${bp[w].n}</span></th>`).join("") + "</tr>";
    for (const p of order) { const tr = document.createElement("tr"); const h = document.createElement("td"); h.className = "rowh"; h.textContent = SHORT[p]; if (NRT.has(p)) h.style.color = "var(--pos)"; tr.appendChild(h);
      for (const w of ws) { const v = bp[w].products[p], td = cell(v.r2, seq(v.r2 / 0.8)); td.title = `90% range ${fmt(v.lo)}–${fmt(v.hi)}`; tr.appendChild(td); }
      t.appendChild(tr); }
    document.getElementById("gauge-period-note").textContent = "r² against the gauges in each 5-year period (hover for the 90% range). Sites with too few values per period are not listed.";
  }
  gs.onchange = drawMonth; ps.onchange = drawPeriod; drawMonth(); drawPeriod();
})();

// ---- gauge check, step by step
(function () {
  const FATE = {used: ["Used", "#1e795f"], fewcells: ["In a cell that never reaches 3 complete gauges", "#e6ab02"],
    nodata: ["No complete month in 2001–2020", "#9aa5a6"], merged: ["Duplicate, merged into a neighbour", "#7570b3"], chirps: ["Used by CHIRPS v3: removed", "#d03b3b"]};
  const cnt = {}; V.stations.forEach(s => cnt[s.fate] = (cnt[s.fate] || 0) + 1);
  const used = V.cells;
  // human labels for the scored cells: site name + letter in map reading order (north to south,
  // west to east within a site), e.g. "Ouémé A". The raw "9.75N 1.75E" stays in hover text only.
  const siteOf = c => c.site.match(/\((.+)\)/)[1];
  const order = [...used.keys()].sort((a, b) => siteOf(used[a]).localeCompare(siteOf(used[b])) || used[b].iy - used[a].iy || used[a].ix - used[b].ix);
  const LABEL = []; { const n = {}; for (const i of order) { const s = siteOf(used[i]); n[s] = (n[s] || 0) + 1; LABEL[i] = `${s} ${String.fromCharCode(64 + n[s])}`; } }
  const ZOOM ={"All": V.bbox, "Niamey (Niger)": [1.4, 12.9, 3.6, 14.6], "Ouémé (Benin)": [1.2, 8.6, 3.2, 10.6], "Gourma (Mali)": [-2.4, 14.4, -0.6, 17.0]};
  let zoom = "All";
  const svg = document.getElementById("gmap");
  function drawGmap() {
    const [x0, y0, x1, y1] = ZOOM[zoom], W = 640, H = Math.round(W * (y1 - y0) / (x1 - x0));
    svg.setAttribute("viewBox", `0 0 ${W} ${H}`);
    const sx = x => (x - x0) / (x1 - x0) * W, sy = y => H - (y - y0) / (y1 - y0) * H;
    const r = zoom === "All" ? 2.6 : 4.5;
    let g = `<rect width="${W}" height="${H}" fill="var(--ocean)"/>`;
    for (const ring of V.countries) g += `<path d="M${ring.xy.map(([x, y]) => sx(x).toFixed(1) + " " + sy(y).toFixed(1)).join("L")}Z" fill="var(--n05)" stroke="var(--n3)" stroke-width="0.8"/>`;
    if (zoom === "All") for (const [n, [x, y]] of Object.entries(V.country_labels)) if (x > x0 && x < x1 && y > y0 && y < y1) g += `<text x="${sx(x)}" y="${sy(y)}" font-size="11" text-anchor="middle" opacity=".6">${n}</text>`;
    used.forEach((c, i) => { const lat = -90 + 0.5 * c.iy, lon = -180 + 0.5 * c.ix;
      g += `<rect class="gcell" data-i="${i}" x="${sx(lon)}" y="${sy(lat + 0.5)}" width="${sx(lon + 0.5) - sx(lon)}" height="${sy(lat) - sy(lat + 0.5)}" fill="var(--b5)" fill-opacity=".12" stroke="var(--b7)" stroke-width="1.5" style="cursor:pointer"><title>${LABEL[i]} (${c.cell}): ${c.months} scored months, median ${c.median_stations_per_month} gauges per month. Click to view in step 5.</title></rect>`;
 });
    for (const f of ["nodata", "merged", "fewcells", "used", "chirps"]) for (const s of V.stations.filter(s => s.fate === f))
      g += `<circle cx="${sx(s.lon).toFixed(1)}" cy="${sy(s.lat).toFixed(1)}" r="${r}" fill="${FATE[f][1]}" stroke="var(--surface)" stroke-width=".6"><title>${s.name} (${s.dataset}), ${s.lat}, ${s.lon}
${s.why}; ${s.complete_months} complete months 2001–2020</title></circle>`;
    // cell letters on top of the dots, with a halo so they stay legible over dense stations
    used.forEach((c, i) => { const lat = -90 + 0.5 * c.iy, lon = -180 + 0.5 * c.ix, big = zoom !== "All";
      g += `<text x="${sx(lon) + (big ? 4 : 2)}" y="${sy(lat + 0.5) + (big ? 14 : 10)}" font-size="${big ? 13 : 10}" font-weight="700" fill="var(--b7)" stroke="var(--surface)" stroke-width="3" paint-order="stroke" pointer-events="none">${LABEL[i].slice(-1)}</text>`; });
    // scale bar ~100 km or 20 km
    const km = zoom === "All" ? 200 : 20, deg = km / (111 * Math.cos((y0 + y1) / 2 * Math.PI / 180));
    g += `<line x1="12" x2="${12 + sx(x0 + deg)}" y1="${H - 14}" y2="${H - 14}" stroke="var(--n9)" stroke-width="2"/><text x="12" y="${H - 20}" font-size="10">${km} km</text>`;
    svg.innerHTML = g;
    svg.querySelectorAll(".gcell").forEach(el => el.addEventListener("click", () => { document.getElementById("ts-cell").value = el.dataset.i; drawTs(); document.getElementById("gts").scrollIntoView({behavior: "smooth", block: "center"}); }));
  }
  document.getElementById("gmap-legend").innerHTML = `<h4>${V.stations.length} gauge records</h4>` + Object.entries(FATE).map(([k, [lab, col]]) =>
    `<div style="display:flex;gap:8px;align-items:center;margin:4px 0"><span class="sw" style="background:${col};border-radius:50%"></span><span style="flex:1">${lab}</span><span class="num">${cnt[k] || 0}</span></div>`).join("") +
    `<div style="display:flex;gap:8px;align-items:center;margin:8px 0 0"><span class="sw" style="background:transparent;border:2px solid var(--b7)"></span><span style="flex:1">Scored 0.5° cell</span><span class="num">${used.length}</span></div>`;
  seg(document.getElementById("gmap-zoom"), Object.keys(ZOOM), zoom, z => { zoom = z; drawGmap(); });
  drawGmap();

  // 0. coverage timeline + dataset table
  (function () {
    const C = V.coverage, ys = C.years, csvg = document.getElementById("coverage"), cw = 17, rh = 18, l = 250, t0 = 24;
    const prodRows = Object.entries(C.products);
    const H = t0 + (C.datasets.length + prodRows.length) * rh + 30, W = l + ys.length * cw + 10;
    csvg.setAttribute("viewBox", `0 0 ${W} ${H}`); csvg.style.minWidth = "760px";
    const mx = Math.max(...C.datasets.flatMap(d => d.per_year));
    let g = "";
    ys.forEach((y, j) => { if (y % 5 === 0) g += `<text x="${l + j * cw + cw / 2}" y="${t0 - 8}" font-size="10" text-anchor="middle">${y}</text>`; });
    g += `<text x="4" y="${t0 + 12}" font-size="10" font-weight="700">GAUGES (complete station-months)</text>`;
    C.datasets.forEach((d, i) => { const y = t0 + (i + 1) * rh;
      g += `<text x="${l - 6}" y="${y + 12}" font-size="11" text-anchor="end">${d.site} (${d.stations})</text>`;
      d.per_year.forEach((v, j) => { if (!v) return; g += `<rect x="${l + j * cw}" y="${y}" width="${cw - 1}" height="${rh - 2}" fill="${rgb(seq(0.2 + 0.8 * Math.sqrt(v / mx)))}"><title>${d.code}, ${ys[j]}: ${v} complete station-months</title></rect>`; }); });
    const yp = t0 + (C.datasets.length + 1) * rh + 8;
    g += `<text x="4" y="${yp + 10}" font-size="10" font-weight="700">PRODUCTS (data available)</text>`;
    prodRows.forEach(([p, [a, b]], i) => { const y = yp + (i + 1) * rh - 6, ya = +a.slice(0, 4), yb = +b.slice(0, 4);
      const j0 = Math.max(0, ys.indexOf(Math.max(ya, ys[0]))), j1 = ys.indexOf(Math.min(yb, ys[ys.length - 1]));
      g += `<text x="${l - 6}" y="${y + 12}" font-size="11" text-anchor="end">${SHORT[p]}</text><rect x="${l + j0 * cw}" y="${y + 4}" width="${(j1 - j0 + 1) * cw - 1}" height="${rh - 8}" fill="${p === "gpcc_full" ? "var(--neg)" : "var(--n3)"}"><title>${SHORT[p]}: ${a} to ${b}</title></rect>`; });
    const w0 = ys.indexOf(C.window[0]), w1 = ys.indexOf(C.window[1]);
    g += `<rect x="${l + w0 * cw - 1}" y="${t0 - 2}" width="${(w1 - w0 + 1) * cw + 1}" height="${H - t0 - 2}" fill="none" stroke="var(--n9)" stroke-width="1.5"/>`;
    g += `<text x="${l + (w1 + 1) * cw + 2}" y="${yp + rh + 6}" font-size="10" fill="var(--neg)">GPCC Full ends 2020</text>`;
    csvg.setAttribute("viewBox", `0 0 ${W + 120} ${H}`);
    csvg.innerHTML = g;
    const tb = document.getElementById("ds-table");
    tb.innerHTML = "<tr><th>Dataset</th><th>Site</th><th>Stations</th><th>Years with data</th><th>Complete station-months 2001–2020</th><th>Stations used</th></tr>" +
      C.datasets.map(d => `<tr><td class="rowh"><code>${d.code}</code></td><td style="text-align:left">${d.site}</td><td class="num">${d.stations}</td><td>${d.first}–${d.last}</td><td class="num">${d.per_year.reduce((s, v, j) => s + (ys[j] >= C.window[0] && ys[j] <= C.window[1] ? v : 0), 0).toLocaleString()}</td><td class="num"><b>${d.used}</b></td></tr>`).join("");
  })();

  // 2. funnel
  const fmax = Math.max(...V.funnel.filter(f => f.unit !== "cell-months").map(f => f.n)), cmax = Math.max(...V.funnel.filter(f => f.unit === "cell-months").map(f => f.n));
  document.getElementById("funnel").innerHTML = V.funnel.map(f => { const cm = f.unit === "cell-months", w = f.n / (cm ? cmax : fmax) * 100;
    return `<span>${f.step}</span><span><span class="fb${cm ? " cm" : ""}" style="width:${w}%">${f.n.toLocaleString()} ${f.unit}</span></span>`; }).join("");

  // 3. independence
  (function () {
    const ov = Object.fromEntries(G.overlap.map(o => [o.cell, o])), rowsC = order.map(i => ov[used[i].cell] && {...ov[used[i].cell], label: LABEL[i]}).filter(Boolean);
    const isvg = document.getElementById("indep"), W = 720, rowH = 22, H = rowsC.length * rowH + 40, l = 150, xmax = Math.max(...rowsC.map(o => o.amma_not_chirps_jjas_median), 4);
    isvg.setAttribute("viewBox", `0 0 ${W} ${H}`);
    const sx = v => l + v / xmax * (W - l - 20);
    let g = "";
    for (let v = 0; v <= xmax; v += 2) g += `<line x1="${sx(v)}" x2="${sx(v)}" y1="18" y2="${H - 18}" stroke="var(--n1)"/><text x="${sx(v)}" y="${H - 4}" font-size="10" text-anchor="middle">${v}</text>`;
    const ser = [["amma_not_chirps_jjas_median", "AMMA gauges used", "#1e795f"], ["gpcc_full_jjas_median", "GPCC Full", "#d03b3b"], ["chirps3_stn_jjas_median", "CHIRPS v3", "#e6ab02"]];
    rowsC.forEach((o, i) => { const y = 24 + i * rowH; g += `<text x="${l - 8}" y="${y + 10}" font-size="11" text-anchor="end">${o.label}</text>`;
      ser.forEach(([k, lab, col], j) => g += `<rect x="${l}" y="${y + j * 6}" width="${Math.max(sx(o[k]) - l, 1)}" height="5" fill="${col}"><title>${lab}: ${o[k]}</title></rect>`); });
    g += ser.map(([k, lab, col], j) => `<rect x="${l + j * 170}" y="2" width="10" height="10" fill="${col}"/><text x="${l + j * 170 + 14}" y="11" font-size="11">${lab}</text>`).join("");
    g += `<text x="${(l + W) / 2}" y="${H - 16}" font-size="10" text-anchor="middle" opacity="0"> </text>`;
    isvg.innerHTML = g;
  })();

  // 4. availability heatmap
  (function () {
    const a = V.availability, asvg = document.getElementById("avail"), cw = 18, ch = 18, l = 110, t = 24, W = l + a.years.length * cw + 10, H = t + used.length * ch + 10;
    asvg.setAttribute("viewBox", `0 0 ${W} ${H}`); asvg.style.minWidth = "620px";
    const mx = Math.max(...a.counts.flat());
    let g = "";
    a.years.forEach((y, j) => { if (y % 5 === 0) g += `<text x="${l + j * cw + cw / 2}" y="${t - 8}" font-size="10" text-anchor="middle">${y}</text>`; });
    order.forEach((i, row) => { const c = used[i]; g += `<text x="${l - 6}" y="${t + row * ch + 13}" font-size="11" text-anchor="end">${LABEL[i]}</text>`;
      a.counts[i].forEach((v, j) => { const col = v ? rgb(seq(0.15 + 0.85 * v / mx)) : "var(--n1)";
        g += `<rect x="${l + j * cw}" y="${t + row * ch}" width="${cw - 1}" height="${ch - 1}" fill="${col}"><title>${LABEL[i]}, ${a.years[j]}: ${v} complete station-months</title></rect>`; }); });
    const j0 = a.years.indexOf(2001), j1 = a.years.indexOf(2020);
    g += `<rect x="${l + j0 * cw - 1}" y="${t - 2}" width="${(j1 - j0 + 1) * cw + 1}" height="${used.length * ch + 3}" fill="none" stroke="var(--n9)" stroke-width="1.5"/>`;
    asvg.innerHTML = g;
  })();

  // 5. time series
  const tsCell = document.getElementById("ts-cell"), tsPer = document.getElementById("ts-period");
  order.forEach(i => { const c = used[i]; tsCell.add(new Option(`${LABEL[i]} (${c.cell}) · ${c.months} scored months`, i)); });
  tsCell.value = used.reduce((b, c, i) => c.months > used[b].months ? i : b, 0);  // start on the cell with most scored months
  [["2001–2005", 0], ["2006–2010", 60], ["2011–2015", 120], ["2016–2020", 180], ["2001–2020", -1]].forEach(([l, v]) => tsPer.add(new Option(l, v)));
  let tsProds = ["chirps_v3", "imerg_late", "era5"];
  const TSCOL = {chirps_v3: "#8a2b1d", imerg_late: "#8e6cc4", era5: "#2a78d6", gpcc_full: "#1b7f5b", imerg_final: "#5a3d99", chirps_v2: "#b5432a"};
  function drawTs() {
    const s = V.series[+tsCell.value], start = +tsPer.value, i0 = start < 0 ? 0 : start, i1 = start < 0 ? V.times.length : start + 60;
    const tsvg = document.getElementById("gts"), W = 760, H = 240, m = {l: 44, r: 10, t: 12, b: 26};
    tsvg.setAttribute("viewBox", `0 0 ${W} ${H}`);
    const idx = []; for (let i = i0; i < i1; i++) idx.push(i);
    const vals = idx.flatMap(i => [s.gauge[i], ...tsProds.map(p => s.products[p][i])]).filter(v => v >= 0);
    const ymax = Math.ceil(Math.max(...vals) / 50) * 50 || 50;
    const sx = i => m.l + (i - i0) / (idx.length - 1) * (W - m.l - m.r), sy = v => H - m.b - v / ymax * (H - m.t - m.b);
    let g = "";
    for (let v = 0; v <= ymax; v += ymax > 300 ? 100 : 50) g += `<line x1="${m.l}" x2="${W - m.r}" y1="${sy(v)}" y2="${sy(v)}" stroke="var(--n1)"/><text x="${m.l - 6}" y="${sy(v) + 3}" font-size="10" text-anchor="end">${v}</text>`;
    idx.forEach(i => { if (V.times[i].endsWith("-01")) g += `<line x1="${sx(i)}" x2="${sx(i)}" y1="${m.t}" y2="${H - m.b}" stroke="var(--n2)"/><text x="${sx(i) + 3}" y="${H - 8}" font-size="10">${V.times[i].slice(0, 4)}</text>`; });
    g += `<text x="10" y="${m.t + 4}" font-size="10" transform="rotate(-90 10 ${(H) / 2}) translate(-${H / 2 - 20} 0)">mm / month</text>`;
    const line = (arr, col, w, dash) => { let d = ""; idx.forEach(i => { const v = arr[i]; if (v < 0) { d += " "; return; } d += (d && !d.endsWith(" ") ? "L" : "M") + sx(i).toFixed(1) + " " + sy(v).toFixed(1); });
      return `<path d="${d.replace(/ +/g, " ").replace(/ M/g, "M")}" fill="none" stroke="${col}" stroke-width="${w}" ${dash ? `stroke-dasharray="${dash}"` : ""}/>`; };
    tsProds.forEach(p => g += line(s.products[p], TSCOL[p], 1.3));
    g += line(s.gauge, "var(--n9)", 2);
    idx.forEach(i => { if (s.ok[i]) g += `<circle cx="${sx(i)}" cy="${sy(s.gauge[i])}" r="2.4" fill="var(--n9)"><title>${V.times[i]}: gauges ${s.gauge[i]} mm (${s.n[i]} stations) · ` + tsProds.map(p => `${SHORT[p]} ${s.products[p][i]}`).join(" · ") + `</title></circle>`; });
    tsvg.innerHTML = g;
  }
  const tsSeg = document.getElementById("ts-prods");
  Object.keys(TSCOL).forEach(p => { const b = document.createElement("button"); b.type = "button"; b.innerHTML = `<span class="sw" style="background:${TSCOL[p]}"></span>${SHORT[p]}`;
    b.setAttribute("aria-pressed", tsProds.includes(p)); b.onclick = () => { tsProds = tsProds.includes(p) ? tsProds.filter(x => x !== p) : [...tsProds, p]; b.setAttribute("aria-pressed", tsProds.includes(p)); drawTs(); }; tsSeg.appendChild(b); });
  tsCell.onchange = drawTs; tsPer.onchange = drawTs; tsPer.value = "180"; drawTs();

  // 6. scatter + 7. dry months
  const scSel = document.getElementById("sc-prod"); V.products.forEach(p => scSel.add(new Option(SHORT[p], p))); scSel.value = "chirps_v3";
  const SITECOL = {"Benin (Ouémé)": "#1e795f", "Niger (Niamey)": "#d95f02", "Mali (Gourma)": "#7570b3"};
  function drawScatter() {
    const p = scSel.value, k = V.products.indexOf(p) + 1, z = V.pooled.z, W = 460, H = 420, m = 40, lim = 3.2;
    const ssvg = document.getElementById("gscatter"); ssvg.setAttribute("viewBox", `0 0 ${W} ${H}`);
    const sx = v => m + (Math.max(-lim, Math.min(lim, v)) + lim) / (2 * lim) * (W - 2 * m), sy = v => H - m - (Math.max(-lim, Math.min(lim, v)) + lim) / (2 * lim) * (H - 2 * m);
    let g = `<rect x="${m}" y="${m}" width="${W - 2 * m}" height="${H - 2 * m}" fill="var(--surface)" stroke="var(--n2)"/>`;
    for (let v = -3; v <= 3; v++) g += `<line x1="${sx(v)}" x2="${sx(v)}" y1="${m}" y2="${H - m}" stroke="var(--n1)"/><line x1="${m}" x2="${W - m}" y1="${sy(v)}" y2="${sy(v)}" stroke="var(--n1)"/><text x="${sx(v)}" y="${H - m + 14}" font-size="10" text-anchor="middle">${v}</text><text x="${m - 6}" y="${sy(v) + 3}" font-size="10" text-anchor="end">${v}</text>`;
    g += `<line x1="${sx(-lim)}" y1="${sy(-lim)}" x2="${sx(lim)}" y2="${sy(lim)}" stroke="var(--n3)" stroke-dasharray="4 3"/>`;
    z.forEach((row, i) => { const site = V.cells[V.pooled.cell[i]].site; g += `<circle cx="${sx(row[0]).toFixed(1)}" cy="${sy(row[k]).toFixed(1)}" r="2" fill="${SITECOL[site]}" fill-opacity=".45"/>`; });
    g += `<text x="${W / 2}" y="${H - 6}" font-size="11" text-anchor="middle">Gauges: how unusual the month was</text><text transform="translate(12 ${H / 2}) rotate(-90)" font-size="11" text-anchor="middle">${SHORT[p]}: how unusual the month was</text>`;
    g += `<text x="${m + 8}" y="${m + 16}" font-size="13" font-weight="700">r² = ${fmt(G.products[p].r2_vs_gauge)}</text>`;
    g += Object.entries(SITECOL).map(([s, c], j) => `<circle cx="${W - m - 110}" cy="${m + 14 + j * 15}" r="4" fill="${c}"/><text x="${W - m - 102}" y="${m + 18 + j * 15}" font-size="10">${s}</text>`).join("");
    ssvg.innerHTML = g;
    // dry-month confusion
    let hit = 0, miss = 0, fa = 0, neither = 0; V.pooled.dry.forEach((d, i) => { if (!V.pooled.valid[i]) return; const gd = d[0], pd = d[k]; if (gd && pd) hit++; else if (gd) miss++; else if (pd) fa++; else neither++; });
    document.getElementById("dry-panel").innerHTML = `<h4>Dry months: ${SHORT[p]} against the gauges</h4><p class="muted" style="margin:0 0 8px">"Dry" = in the driest third of years for that cell and calendar month, each source judged on its own record.</p>
      <table class="t" style="font-size:12px"><tr><th></th><th>Gauges dry</th><th>Gauges not dry</th></tr>
      <tr><td class="rowh">${SHORT[p]} dry</td><td style="background:var(--b05)"><b>${hit}</b> hits</td><td>${fa} false alarms</td></tr>
      <tr><td class="rowh">${SHORT[p]} not dry</td><td>${miss} misses</td><td>${neither}</td></tr></table>
      <p style="margin:8px 0 0">Hit rate ${pct(hit / (hit + miss))} (chance 33%) · false alarms ${pct(fa / (hit + fa))} of its dry calls.</p>`;
  }
  scSel.onchange = drawScatter; drawScatter();
})();

// ---- gauge 3-month and TAHMO tables
(function () {
  const ci = c => `<span class="muted num" style="font-size:11px"> ${fmt(c[0])}–${fmt(c[1])}</span>`;
  function table(el, grp, prods) {
    const ord = [...prods].sort((a, b) => (grp.products[b].r2 ?? -1) - (grp.products[a].r2 ?? -1));
    el.innerHTML = "<tr><th>Product</th><th>r²</th><th>Spearman</th><th>Dry hit rate</th><th>False alarms</th></tr>";
    for (const p of ord) { const v = grp.products[p], tr = document.createElement("tr");
      const h = document.createElement("td"); h.className = "rowh"; h.textContent = SHORT[p]; if (NRT.has(p)) h.style.color = "var(--pos)"; tr.appendChild(h);
      const r = cell(v.r2, seq(v.r2 / 0.7)); r.innerHTML += ci(v.r2_ci); tr.appendChild(r);
      tr.appendChild(cell(v.spearman, seq(v.spearman / 0.85)));
      tr.appendChild(cell(v.dry_hit, seq((v.dry_hit - 0.33) / 0.45)));
      const fa = document.createElement("td"); fa.className = "num"; fa.textContent = v.dry_false_alarm == null ? "–" : pct(v.dry_false_alarm); tr.appendChild(fa);
      el.appendChild(tr); }
  }
  const s3 = document.getElementById("g3-site");
  Object.keys(G3.groups).forEach(k => s3.add(new Option(k === "All" ? `All gauged cells (${G3.groups[k].n} windows)` : `${k} (${G3.groups[k].n} windows)`, k)));
  const draw3 = () => table(document.getElementById("gauge3-table"), G3.groups[s3.value], G3.products);
  s3.onchange = draw3; draw3();
  let tScale = "3-month";
  function drawT() {
    const src = tScale === "monthly" ? T.groups : T.groups_3month, el = document.getElementById("tahmo-table");
    const regs = ["All|All", "East Africa|All", "West Africa|All", "Southern Africa|All", "Central Africa|All"].filter(k => src[k]);
    const prods = T.products, order = [...prods].sort((a, b) => (src["All|All"].products[b].r2 ?? -1) - (src["All|All"].products[a].r2 ?? -1));
    el.innerHTML = "<tr><th>Product</th>" + regs.map(k => `<th>${k.split("|")[0] === "All" ? "All regions" : k.split("|")[0]}<br><span class="muted" style="font-weight:400">${src[k].cells} cells · r² / dry hit</span></th>`).join("") + "</tr>";
    for (const p of order) { const tr = document.createElement("tr"); const h = document.createElement("td"); h.className = "rowh"; h.textContent = SHORT[p]; if (NRT.has(p)) h.style.color = "var(--pos)"; tr.appendChild(h);
      for (const k of regs) { const v = src[k].products[p], td = cell(v.r2, seq(v.r2 / 0.5)); td.innerHTML += `<span class="num" style="font-size:11px;opacity:.85"> / ${v.dry_hit == null ? "–" : pct(v.dry_hit)}</span>`;
        td.title = `r² 90% range ${fmt(v.r2_ci[0])}–${fmt(v.r2_ci[1])}; Spearman ${fmt(v.spearman)}`; tr.appendChild(td); }
      el.appendChild(tr); }
  }
  seg(document.getElementById("t-scale"), ["3-month", "monthly"], tScale, s => { tScale = s; drawT(); }); drawT();
})();

// ---- TAHMO selection funnel (counts only: no station locations or per-country records on the public page)
(function () {
  const fmax = Math.max(...TV.funnel.filter(f => f.unit === "stations").map(f => f.n)), cmax = Math.max(...TV.funnel.filter(f => f.unit !== "stations").map(f => f.n));
  document.getElementById("tfunnel").innerHTML = TV.funnel.map(f => { const cm = f.unit !== "stations", w = f.n / (cm ? cmax : fmax) * 100;
    return `<span>${f.step}</span><span><span class="fb${cm ? " cm" : ""}" style="width:${w}%">${f.n.toLocaleString()} ${f.unit}</span></span>`; }).join("");
  document.getElementById("tfunnel").id = "tfunnel"; document.getElementById("tfunnel").style.cssText = "display:grid;grid-template-columns:auto minmax(0,1fr);gap:5px 12px;align-items:center;font-size:13px;margin:8px 0";
})();

// ---- initial draw (main result: CHIRP v3 (no stations) and IMERG Late share mistakes)
drawInputs(); drawDep(); drawSkill(); drawDecision();
const redraw = () => { drawDep(); drawSkill(); drawDecision(); drawMap(); };
matchMedia("(prefers-color-scheme: dark)").addEventListener?.("change", redraw);
new MutationObserver(redraw).observe(document.documentElement, {attributes: true, attributeFilter: ["data-theme"]});
return initMap;
}
