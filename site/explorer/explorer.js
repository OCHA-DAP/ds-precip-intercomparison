// Rainfall product correlation explorer. Everything is computed in the browser
// from per-pixel yearly totals (data/<product>_<agg>.bin, see scripts/explorer_data.py).

const DIV = ["#8f1d1d", "#d03b3b", "#ef9a8a", "#f0efec", "#86b6ef", "#2a78d6", "#104281"];
const NODATA = [215, 219, 219, 255];
const LAT0 = -60, LAT1 = 84, RES = 0.5;
const $ = (s) => document.querySelector(s);
const fmt = (v, d = 2) => (v == null || !isFinite(v) ? "–" : (+v).toFixed(d));

let META, CELLS, GRID, UNITS = { country: new Map(), admin1: new Map() };
const cache = new Map();
const state = { a: "chirps_v2", b: "era5", agg: "annual", stat: "pearson", unit: "pixel", cell: null, country: null, admin1: null };
let map, overlay, corrField = null, nField = null, sel = { A: null, B: null, years: [] };

// ------------------------------------------------------------------ data
async function getBin(path) {
  const r = await fetch(path);
  if (!r.ok) throw new Error(`${path}: HTTP ${r.status}`);
  return r.arrayBuffer();
}

async function loadCells() {
  const buf = await getBin("data/cells.bin");
  const N = META.ncell;
  let off = 0;
  const take = (Type, n) => { const a = new Type(buf.slice(off, off + n * Type.BYTES_PER_ELEMENT)); off += n * Type.BYTES_PER_ELEMENT; return a; };
  CELLS = {
    iy: take(Int16Array, N), ix: take(Int16Array, N), country: take(Uint16Array, N), admin1: take(Uint16Array, N),
    w: take(Float32Array, N), rainy: take(Uint16Array, N), start: take(Uint8Array, N),
  };
  GRID = new Int32Array(720 * 360).fill(-1);
  for (let k = 0; k < N; k++) GRID[CELLS.iy[k] * 720 + CELLS.ix[k]] = k;
  for (let k = 0; k < N; k++) {
    for (const u of ["country", "admin1"]) {
      const v = CELLS[u][k];
      if (v === 65535) continue;
      if (!UNITS[u].has(v)) UNITS[u].set(v, []);
      UNITS[u].get(v).push(k);
    }
  }
}

async function series(product, agg) {
  const key = `${product}_${agg}`;
  if (cache.has(key)) return cache.get(key);
  const p = (async () => {
    const f = META.files[key];
    const buf = await getBin(`data/${key}.bin`);
    const N = META.ncell, Y = f.nyear;
    const q = new Uint8Array(buf, 0, N * Y);
    const lo = new Uint16Array(buf.slice(N * Y, N * Y + 2 * N));
    const hi = new Uint16Array(buf.slice(N * Y + 2 * N, N * Y + 4 * N));
    const v = new Float32Array(N * Y);
    for (let k = 0; k < N; k++) {
      const a = lo[k], s = (hi[k] - lo[k]) / 254;
      for (let y = 0; y < Y; y++) { const x = q[k * Y + y]; v[k * Y + y] = x === 255 ? NaN : a + x * s; }
    }
    return { first: f.first_year, n: Y, v };
  })();
  cache.set(key, p);
  return p;
}

// ------------------------------------------------------------------ statistics
function ranks(arr) {
  const idx = arr.map((v, i) => i).sort((i, j) => arr[i] - arr[j]);
  const r = new Array(arr.length);
  for (let i = 0; i < idx.length;) {
    let j = i; while (j + 1 < idx.length && arr[idx[j + 1]] === arr[idx[i]]) j++;
    for (let k = i; k <= j; k++) r[idx[k]] = (i + j) / 2 + 1;
    i = j + 1;
  }
  return r;
}
function detrend(t, x) {
  const n = x.length, mt = t.reduce((a, b) => a + b, 0) / n, mx = x.reduce((a, b) => a + b, 0) / n;
  let sxy = 0, sxx = 0;
  for (let i = 0; i < n; i++) { sxy += (t[i] - mt) * (x[i] - mx); sxx += (t[i] - mt) ** 2; }
  const b = sxx > 0 ? sxy / sxx : 0;
  return x.map((v, i) => v - mx - b * (t[i] - mt));
}
function pearson(x, y) {
  const n = x.length; if (n < 3) return NaN;
  let mx = 0, my = 0; for (let i = 0; i < n; i++) { mx += x[i]; my += y[i]; } mx /= n; my /= n;
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) { const dx = x[i] - mx, dy = y[i] - my; sxy += dx * dy; sxx += dx * dx; syy += dy * dy; }
  return sxx > 0 && syy > 0 ? sxy / Math.sqrt(sxx * syy) : NaN;
}
function corr(t, x, y, stat) {
  if (stat === "spearman") return pearson(ranks(x), ranks(y));
  if (stat === "detrended") return pearson(detrend(t, x), detrend(t, y));
  return pearson(x, y);
}
// paired values of cell k over the overlapping years
function pairAt(A, B, k, y0, n) {
  const t = [], x = [], y = [];
  for (let i = 0; i < n; i++) {
    const yr = y0 + i, a = A.v[k * A.n + (yr - A.first)], b = B.v[k * B.n + (yr - B.first)];
    if (isFinite(a) && isFinite(b)) { t.push(yr); x.push(a); y.push(b); }
  }
  return { t, x, y };
}
function overlap(A, B) {
  const y0 = Math.max(A.first, B.first), y1 = Math.min(A.first + A.n, B.first + B.n);
  return { y0, n: Math.max(0, y1 - y0) };
}

// ------------------------------------------------------------------ map
function colorFor(r) {
  if (!isFinite(r)) return null;
  const edges = [-0.5, 0, 0.2, 0.4, 0.6, 0.8]; // 7 classes
  let i = 0; while (i < edges.length && r > edges[i]) i++;
  return DIV[i];
}
const hex = (h) => [parseInt(h.slice(1, 3), 16), parseInt(h.slice(3, 5), 16), parseInt(h.slice(5, 7), 16), 255];
const mercY = (lat) => Math.log(Math.tan(Math.PI / 4 + (lat * Math.PI) / 360));

function render() {
  const W = 1440; // 2 px per 0.5 deg cell
  const y0 = mercY(LAT0), y1 = mercY(LAT1);
  const H = Math.round(((y1 - y0) / (2 * Math.PI)) * W);
  const cv = document.createElement("canvas"); cv.width = W; cv.height = H;
  const ctx = cv.getContext("2d"); const img = ctx.createImageData(W, H);
  const rgba = DIV.map(hex);
  for (let j = 0; j < H; j++) {
    const my = y1 - ((j + 0.5) / H) * (y1 - y0);
    const lat = (2 * Math.atan(Math.exp(my)) - Math.PI / 2) * 180 / Math.PI;
    const iy = Math.floor((lat + 90) / RES);
    for (let i = 0; i < W; i++) {
      const ix = Math.floor(i / 2), k = GRID[iy * 720 + ix];
      if (k < 0) continue;
      const r = corrField[k], o = (j * W + i) * 4;
      let c;
      if (isFinite(r) && nField[k] >= 10) {
        const col = colorFor(r); c = rgba[DIV.indexOf(col)];
      } else c = NODATA;
      img.data[o] = c[0]; img.data[o + 1] = c[1]; img.data[o + 2] = c[2]; img.data[o + 3] = 235;
    }
  }
  ctx.putImageData(img, 0, 0);
  const url = cv.toDataURL();
  if (overlay) overlay.setUrl(url);
  else overlay = L.imageOverlay(url, [[LAT0, -180], [LAT1, 180]], { interactive: false, pane: "overlayPane" }).addTo(map);
}

async function update() {
  $("#status").textContent = "loading…";
  try {
    const [A, B] = await Promise.all([series(state.a, state.agg), series(state.b, state.agg)]);
    sel = { A, B };
    const { y0, n } = overlap(A, B);
    const N = META.ncell;
    corrField = new Float32Array(N); nField = new Uint8Array(N);
    const vals = [];
    for (let k = 0; k < N; k++) {
      const p = pairAt(A, B, k, y0, n);
      nField[k] = p.t.length;
      corrField[k] = p.t.length >= 10 ? corr(p.t, p.x, p.y, state.stat) : NaN;
      if (isFinite(corrField[k])) vals.push(corrField[k]);
    }
    vals.sort((a, b) => a - b);
    const med = vals.length ? vals[Math.floor(vals.length / 2)] : NaN;
    $("#map-summary").textContent = n ? `median over ${vals.length.toLocaleString()} pixels: ${fmt(med)} · years ${y0}–${y0 + n - 1}` : "no overlapping years";
    render();
    drawScatter();
    $("#status").textContent = "";
    writeHash();
  } catch (e) {
    console.error(e); $("#status").textContent = `failed: ${e.message}`;
  }
}

// ------------------------------------------------------------------ scatter
function unitSeries() {
  const { A, B } = sel; if (!A || !B) return null;
  const { y0, n } = overlap(A, B);
  if (state.unit === "pixel") {
    if (state.cell == null) return null;
    const p = pairAt(A, B, state.cell, y0, n);
    const lat = CELLS.iy[state.cell] * RES - 90 + RES / 2, lon = CELLS.ix[state.cell] * RES - 180 + RES / 2;
    const ci = CELLS.country[state.cell], ai = CELLS.admin1[state.cell];
    const where = [ai !== 65535 ? META.admin1[ai].name : null, ci !== 65535 ? META.countries[ci].name : null].filter(Boolean).join(", ");
    return { ...p, title: `Pixel ${fmt(lat, 2)}°, ${fmt(lon, 2)}°${where ? ` — ${where}` : ""}`, cells: [state.cell] };
  }
  const id = state.unit === "country" ? state.country : state.admin1;
  if (id == null) return null;
  const cells = UNITS[state.unit].get(id) || [];
  const t = [], x = [], y = [];
  for (let i = 0; i < n; i++) {
    const yr = y0 + i; let sa = 0, sb = 0, sw = 0, wall = 0;
    for (const k of cells) {
      const w = CELLS.w[k]; wall += w;
      const a = A.v[k * A.n + (yr - A.first)], b = B.v[k * B.n + (yr - B.first)];
      if (isFinite(a) && isFinite(b)) { sa += w * a; sb += w * b; sw += w; }
    }
    if (sw > 0 && sw >= 0.5 * wall) { t.push(yr); x.push(sa / sw); y.push(sb / sw); }
  }
  const name = state.unit === "country" ? META.countries[id].name : `${META.admin1[id].name}, ${META.admin1[id].iso3}`;
  return { t, x, y, title: `${name} (${cells.length} cell${cells.length === 1 ? "" : "s"})`, cells };
}

function niceTicks(lo, hi, n) {
  const span = hi - lo || 1, step0 = span / n, mag = Math.pow(10, Math.floor(Math.log10(step0)));
  const step = [1, 2, 2.5, 5, 10].map((k) => k * mag).find((s) => span / s <= n) || 10 * mag;
  const out = []; for (let v = Math.ceil(lo / step) * step; v <= hi + 1e-9; v += step) out.push(+v.toFixed(6));
  return out;
}
const svgEl = (tag, attrs = {}) => { const e = document.createElementNS("http://www.w3.org/2000/svg", tag); for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v); return e; };
const label = (k) => META.products.find((p) => p.key === k)?.label.replace(/ \(.*\)$/, "") || k;
const aggLabel = () => META.aggs.find((a) => a.key === state.agg)?.label || state.agg;

function drawScatter() {
  const host = $("#scatter"); host.innerHTML = "";
  const s = unitSeries();
  if (!s || s.t.length === 0) {
    $("#scatter-title").textContent = s ? `${s.title}: no overlapping years` : (state.unit === "pixel" ? "Click the map to pick a pixel" : "Pick a unit above");
    $("#scatter-stats").innerHTML = ""; return;
  }
  $("#scatter-title").textContent = s.title;
  const W = 460, H = 420, m = { l: 58, r: 12, t: 10, b: 46 };
  const all = s.x.concat(s.y); let lo = Math.min(...all), hi = Math.max(...all);
  const pad = (hi - lo) * 0.08 || 10; lo = Math.max(0, lo - pad); hi += pad;
  const X = (v) => m.l + ((v - lo) / (hi - lo)) * (W - m.l - m.r);
  const Y = (v) => H - m.b - ((v - lo) / (hi - lo)) * (H - m.t - m.b);
  const svg = svgEl("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": `${label(state.a)} vs ${label(state.b)}` });
  const g = svgEl("g", { class: "grid" }), ax = svgEl("g", { class: "axis" });
  for (const tv of niceTicks(lo, hi, 5)) {
    g.append(svgEl("line", { x1: X(tv), x2: X(tv), y1: m.t, y2: H - m.b }), svgEl("line", { x1: m.l, x2: W - m.r, y1: Y(tv), y2: Y(tv) }));
    const tx = svgEl("text", { x: X(tv), y: H - m.b + 15, "text-anchor": "middle" }); tx.textContent = tv; ax.append(tx);
    const ty = svgEl("text", { x: m.l - 6, y: Y(tv) + 4, "text-anchor": "end" }); ty.textContent = tv; ax.append(ty);
  }
  const xl = svgEl("text", { x: m.l + (W - m.l - m.r) / 2, y: H - 8, "text-anchor": "middle" }); xl.textContent = `${label(state.a)} (mm)`; ax.append(xl);
  const yl = svgEl("text", { x: 14, y: m.t + (H - m.t - m.b) / 2, transform: `rotate(-90 14 ${m.t + (H - m.t - m.b) / 2})`, "text-anchor": "middle" }); yl.textContent = `${label(state.b)} (mm)`; ax.append(yl);
  svg.append(g, ax);
  svg.append(svgEl("line", { x1: X(lo), y1: Y(lo), x2: X(hi), y2: Y(hi), stroke: "#5e6a6b", "stroke-dasharray": "4 3", "stroke-width": 1 }));
  const t0 = Math.min(...s.t), t1 = Math.max(...s.t);
  s.t.forEach((yr, i) => {
    const f = t1 > t0 ? (yr - t0) / (t1 - t0) : 1;
    const col = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281"][Math.min(4, Math.floor(f * 5))];
    const tx = svgEl("text", { x: X(s.x[i]), y: Y(s.y[i]) + 3.5, "text-anchor": "middle", class: "yr", fill: col });
    tx.textContent = `'${String(yr).slice(2)}`;
    const tt = svgEl("title"); tt.textContent = `${yr}: ${label(state.a)} ${fmt(s.x[i], 0)} mm, ${label(state.b)} ${fmt(s.y[i], 0)} mm`;
    tx.append(tt); svg.append(tx);
  });
  host.append(svg);
  const r = pearson(s.x, s.y), rho = pearson(ranks(s.x), ranks(s.y)), rd = pearson(detrend(s.t, s.x), detrend(s.t, s.y));
  const ratio = s.y.reduce((a, b) => a + b, 0) / s.x.reduce((a, b) => a + b, 0);
  let extra = "";
  if (state.unit === "pixel") {
    const mask = CELLS.rainy[state.cell], st = CELLS.start[state.cell], M = "JFMAMJJASOND";
    const months = []; for (let i = 0; i < 12; i++) { const mo = (st + i) % 12; if (mask >> mo & 1) months.push(M[mo]); }
    extra = `<br>Rainy-season months here: <b>${months.join(" ")}</b> (hydrological year starts in ${"Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec".split(" ")[st]})`;
  }
  $("#scatter-stats").innerHTML = `${aggLabel()} · ${s.t.length} years (${t0}–${t1}) · Pearson r <b>${fmt(r)}</b> · Spearman ρ ${fmt(rho)} ·
    detrended r ${fmt(rd)} · ${label(state.b)} ÷ ${label(state.a)} = ${fmt(ratio)}<br>Labels are years ('95 = 1995), darker = more recent;
    dashed line = 1:1.${extra}`;
}

// ------------------------------------------------------------------ UI
function writeHash() {
  const p = new URLSearchParams({ a: state.a, b: state.b, agg: state.agg, stat: state.stat, unit: state.unit });
  if (state.cell != null) { p.set("lat", fmt(CELLS.iy[state.cell] * RES - 90 + RES / 2, 2)); p.set("lon", fmt(CELLS.ix[state.cell] * RES - 180 + RES / 2, 2)); }
  if (state.unit === "country" && state.country != null) p.set("c", META.countries[state.country].iso3);
  if (state.unit === "admin1" && state.admin1 != null) p.set("a1", String(state.admin1));
  history.replaceState(null, "", `#${p}`);
}
let HASH = new URLSearchParams();
function readHash() {
  const p = new URLSearchParams(location.hash.slice(1));
  HASH = p;
  const valid = {
    a: META.products.map((x) => x.key), b: META.products.map((x) => x.key), agg: META.aggs.map((x) => x.key),
    stat: ["pearson", "spearman", "detrended"], unit: ["pixel", "admin1", "country"],
  };
  for (const k of Object.keys(valid)) if (valid[k].includes(p.get(k))) state[k] = p.get(k);
  if (p.get("a1") != null && META.admin1[+p.get("a1")]) state.admin1 = +p.get("a1");
  if (p.get("lat") && p.get("lon")) state.cell = cellAt(+p.get("lat"), +p.get("lon"));
  if (p.get("c")) { const i = META.countries.findIndex((c) => c.iso3 === p.get("c")); if (i >= 0) state.country = i; }
  if (state.admin1 != null && state.country == null) {
    const i = META.countries.findIndex((c) => c.iso3 === META.admin1[state.admin1].iso3); if (i >= 0) state.country = i;
  }
}
function cellAt(lat, lon) {
  const iy = Math.floor((lat + 90) / RES), ix = Math.floor((((lon + 180) % 360) + 360) % 360 / RES);
  if (iy < 0 || iy >= 360) return null;
  const k = GRID[iy * 720 + ix]; return k >= 0 ? k : null;
}
function pickCell(k, setUnits = true) {
  if (k == null) return;
  state.cell = k;
  if (setUnits && CELLS.country[k] !== 65535) { state.country = CELLS.country[k]; $("#sel-country").value = String(state.country); }
  if (setUnits && CELLS.admin1[k] !== 65535) { state.admin1 = CELLS.admin1[k]; fillAdmin1(); $("#sel-admin1").value = String(state.admin1); }
  if (marker) marker.setLatLng([CELLS.iy[k] * RES - 90 + RES / 2, CELLS.ix[k] * RES - 180 + RES / 2]);
  else marker = L.circleMarker([CELLS.iy[k] * RES - 90 + RES / 2, CELLS.ix[k] * RES - 180 + RES / 2], { radius: 6, color: "#1f2324", weight: 2, fill: false }).addTo(map);
  drawScatter(); writeHash();
}
let marker = null;
function fillAdmin1() {
  const s = $("#sel-admin1"); s.innerHTML = "";
  const iso = state.country != null ? META.countries[state.country].iso3 : null;
  const opts = META.admin1.map((a, i) => ({ a, i })).filter(({ a, i }) => (!iso || a.iso3 === iso) && UNITS.admin1.has(i));
  s.append(new Option(iso ? "Admin-1 unit…" : "Admin-1 unit (pick a country first)…", ""));
  opts.sort((p, q) => p.a.name.localeCompare(q.a.name)).forEach(({ a, i }) => s.append(new Option(a.name, String(i))));
}

(async function main() {
  try {
    META = await (await fetch("data/meta.json")).json();
    await loadCells();
    readHash();
    for (const id of ["#sel-a", "#sel-b"]) META.products.forEach((p) => $(id).append(new Option(p.label, p.key)));
    META.aggs.forEach((a) => $("#sel-agg").append(new Option(a.label, a.key)));
    $("#sel-a").value = state.a; $("#sel-b").value = state.b; $("#sel-agg").value = state.agg; $("#sel-stat").value = state.stat;
    const cs = $("#sel-country"); cs.append(new Option("Country…", ""));
    META.countries.map((c, i) => ({ c, i })).filter(({ i }) => UNITS.country.has(i)).sort((p, q) => p.c.name.localeCompare(q.c.name))
      .forEach(({ c, i }) => cs.append(new Option(c.name, String(i))));
    if (state.country != null) cs.value = String(state.country);
    fillAdmin1();
    if (state.admin1 != null) $("#sel-admin1").value = String(state.admin1);
    const seg = $("#seg-unit");
    [["pixel", "Pixel"], ["admin1", "Admin 1"], ["country", "Country"]].forEach(([v, lab]) => {
      const b = document.createElement("button"); b.type = "button"; b.textContent = lab; b.setAttribute("aria-pressed", String(state.unit === v));
      b.addEventListener("click", () => { state.unit = v; seg.querySelectorAll("button").forEach((x) => x.setAttribute("aria-pressed", String(x === b))); drawScatter(); writeHash(); });
      seg.append(b);
    });
    $("#sel-a").addEventListener("change", (e) => { state.a = e.target.value; update(); });
    $("#sel-b").addEventListener("change", (e) => { state.b = e.target.value; update(); });
    $("#sel-agg").addEventListener("change", (e) => { state.agg = e.target.value; update(); });
    $("#sel-stat").addEventListener("change", (e) => { state.stat = e.target.value; update(); });
    cs.addEventListener("change", () => {
      state.country = cs.value === "" ? null : +cs.value;
      state.admin1 = null; fillAdmin1();
      state.unit = "country";
      document.querySelectorAll("#seg-unit button").forEach((x) => x.setAttribute("aria-pressed", String(x.textContent === "Country")));
      drawScatter(); writeHash();
    });
    $("#sel-admin1").addEventListener("change", (e) => { state.admin1 = e.target.value === "" ? null : +e.target.value; if (state.admin1 != null) { state.unit = "admin1"; document.querySelectorAll("#seg-unit button").forEach((x) => x.setAttribute("aria-pressed", String(x.textContent === "Admin 1"))); } drawScatter(); writeHash(); });
    $("#legend-bar").style.background = `linear-gradient(to right, ${DIV.map((c, i) => `${c} ${(i / DIV.length) * 100}%, ${c} ${((i + 1) / DIV.length) * 100}%`).join(", ")})`;

    map = L.map("map", { worldCopyJump: true, minZoom: 1, maxZoom: 10 }).setView([10, 15], 2);
    map.attributionControl.setPrefix('<a href="https://leafletjs.com">Leaflet</a>');
    L.tileLayer("https://server.arcgisonline.com/ArcGIS/rest/services/Canvas/World_Light_Gray_Base/MapServer/tile/{z}/{y}/{x}",
      { attribution: "Tiles © Esri — Esri, HERE, Garmin, OpenStreetMap contributors", maxZoom: 16 }).addTo(map);
    map.createPane("labels"); map.getPane("labels").style.zIndex = 650; map.getPane("labels").style.pointerEvents = "none";
    L.tileLayer("https://server.arcgisonline.com/ArcGIS/rest/services/Canvas/World_Light_Gray_Reference/MapServer/tile/{z}/{y}/{x}",
      { pane: "labels", maxZoom: 16 }).addTo(map);
    const tip = L.tooltip({ className: "ex-tip", direction: "top", offset: [0, -6] });
    map.on("mousemove", (e) => {
      const k = cellAt(e.latlng.lat, e.latlng.lng);
      if (k == null || !corrField) { map.closeTooltip(tip); return; }
      const why = nField[k] < 10 ? "too few overlapping years" : "no year-to-year variation (e.g. a rainless season)";
      tip.setLatLng(e.latlng).setContent(isFinite(corrField[k]) ? `r = ${fmt(corrField[k])} (${nField[k]} years)` : `no value: ${why} (${nField[k]} years)`);
      map.openTooltip(tip);
    });
    map.on("mouseout", () => map.closeTooltip(tip));
    map.on("click", (e) => { const k = cellAt(e.latlng.lat, e.latlng.lng); if (k != null) { if (state.unit !== "pixel") { state.unit = "pixel"; document.querySelectorAll("#seg-unit button").forEach((x) => x.setAttribute("aria-pressed", String(x.textContent === "Pixel"))); } pickCell(k); } });
    // start from the linked pixel/unit; otherwise a Niamey pixel
    const linkedUnit = HASH.get("c") != null || HASH.get("a1") != null;
    if (state.cell == null && !linkedUnit) state.cell = cellAt(13.25, 2.25);
    await update();
    if (state.cell != null) pickCell(state.cell, !linkedUnit);
  } catch (e) {
    console.error(e); const f = $("#fatal"); f.hidden = false; f.textContent = `Could not load the explorer data (${e.message}).`;
  }
})();
