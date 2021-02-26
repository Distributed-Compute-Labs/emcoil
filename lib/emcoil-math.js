 function workerfn(...args) {
 (function(){function r(e,n,t){function o(i,f){if(!n[i]){if(!e[i]){var c="function"==typeof require&&require;if(!f&&c)return c(i,!0);if(u)return u(i,!0);var a=new Error("Cannot find module '"+i+"'");throw a.code="MODULE_NOT_FOUND",a}var p=n[i]={exports:{}};e[i][0].call(p.exports,function(r){var n=e[i][1][r];return o(n||r)},p,p.exports,r,e,n,t)}return n[i].exports}for(var u="function"==typeof require&&require,i=0;i<t.length;i++)o(t[i]);return o}return r})()({1:[function(require,module,exports){

// Import from SparcMath
let SparcMath = require('../sparc-math/_es2015/sparc-math.js');
// console.log(SparcMath);
let besselI = SparcMath.default.besseli;
let struveL = SparcMath.default.struvel;
let Complex = SparcMath.default.Complex;
let quad = SparcMath.default.quad;
let amos = SparcMath.default.amos;
let besselK = amos.besselK;

if (typeof window !== 'undefined') {
 window.SparcMath = SparcMath.default
}

if (typeof window === 'undefined' && typeof self !== 'undefined'){
 window = self
}

const _$_mu0 = 0;
const _$_n = 0;
const _$_i = 0;
const _$_wc = 0;
const _$_wa = 0;
const _$_R = 0;
const _$_mu1 = 0;
const _$_sigma1 = 0;
const _$_w1 = 0;
const _$_mu2 = 0;
const _$_sigma2 = 0;
const _$_w2 = 0;
const _$_f = 0;
const _$_v = 0;

const mu0 = _$_mu0 || 1.25663706e-6 // - Vacuum permeability
// const n = _$_n || 9.46e6 // m^-2 - Coil wire turn density
// const i = _$_i || 1; // A - Coil current
// const wc = _$_wc || 0.005; // m - Coil width
// const wa = _$_wa || 0.005; // m - Air gap between coil and OUTER tube
// const R = _$_R || 0.05; // m - Outer radius of OUTER tube
// const mu1 = _$_mu1 || 200; // - Relative magnetic permeability of OUTER tube
// const sigma1 = _$_sigma1 || 1.04e7; // S * m^-1 - Conductivity of OUTER tube
// const w1 = _$_w1 || 0.005; // m - Thickness of OUTER tube
// const mu2 = _$_mu2 || 0.9996 // - Relative magnetic permeability of INNER tube
// const sigma2 = _$_sigma2 || 2.22; // S * m^-1 - Conductivity of INNER tube
// const w2 = _$_w2 || 0.005; // m - Thickness of INNER tube
// const f = _$_f || 1000; // Hz - Frequency of alternating current
// const v = _$_v || 10; // m * s^-1 - Coil velocity


// Coil function sub-functions
// function length(R, wa, wc) {
// return (0.000013 / (wc * (2 * (R + wa) + wc)))
// }

function struveL_besselK(k, ri, ro) {
 let term1 = Complex.mul(besselK(1, k * ro), struveL(0, k * ro) * ro);
 let term2 = Complex.mul(besselK(0, k * ro), struveL(1, k * ro) * ro);
 let term3 = Complex.mul(besselK(1, k * ri), struveL(0, k * ri) * ri);
 let term4 = Complex.mul(besselK(0, k * ri), struveL(1, k * ri) * ri);
 return Complex.sub(Complex.add(term1, term2), Complex.add(term3, term4));
}

//Tube function sub-functions
function Lambda1(k, sign, mu1, sigma1, f, v) {
 return Complex.sqrt([k**2, (2 * Math.PI * f + sign * k * v) * mu0 * mu1 * sigma1]);
}

function Lambda2(k, sign, mu2, sigma2, f, v) {
 return Complex.sqrt([k**2, (2 * Math.PI * f + sign * k * v) * mu0 * mu2 * sigma2]);
}

function scriptA(k, Lambda1, R, mu1) {
 let Lambda1_a = Complex.mul(Lambda1, R);
 let term1 = Complex.prod([k, besselK(0, k * R), besselI(1, Lambda1_a)]);
 let term2 = Complex.prod([mu1, Lambda1, besselI(0, Lambda1_a), besselK(1, k * R)]);
 return Complex.add(term1, term2);
}

function scriptB(k, Lambda1, R, mu1) {
 let Lambda1_a = Complex.mul(Lambda1, R);
 let term1 = Complex.prod([k, besselK(0, k * R), besselK(1, Lambda1_a)]);
 let term2 = Complex.prod([mu1, Lambda1, besselK(0, Lambda1_a), besselK(1, k * R)]);
 return Complex.sub(term1, term2);
}

function scriptC(k, Lambda1, R, mu1) {
 let Lambda1_a = Complex.mul(Lambda1, R);
 let term1 = Complex.prod([k, besselI(0, k * R), besselI(1, Lambda1_a)]);
 let term2 = Complex.prod([mu1, Lambda1, besselI(0, Lambda1_a), besselI(1, k * R)]);
 return Complex.sub(term1, term2);
}

function scriptD(k, Lambda1, R, mu1) {
 let a = R
 let Lambda1_a = Complex.mul(Lambda1, a);
 let term1 = Complex.prod([k, besselI(0, k * a), besselK(1, Lambda1_a)]);
 let term2 = Complex.prod([mu1, Lambda1, besselK(0, Lambda1_a), besselI(1, k * a)]);
 return Complex.add(term1, term2);
}

function scriptE(k, Lambda1, Lambda2, R, w1, mu1, mu2) {
 let b = R - w1;
 let Lambda1_b = Complex.mul(Lambda1, b);
 let Lambda2_b = Complex.mul(Lambda2, b);
 let term1 = Complex.prod([mu1, Lambda1, besselK(0, Lambda1_b), besselI(1, Lambda2_b)]);
 let term2 = Complex.prod([mu2, Lambda2, besselI(0, Lambda2_b), besselK(1, Lambda1_b)]);
 return Complex.add(term1, term2);
}

function scriptF(k, Lambda1, Lambda2, R, w1, mu1, mu2) {
 let b = R - w1;
 let Lambda1_b = Complex.mul(Lambda1, b);
 let Lambda2_b = Complex.mul(Lambda2, b);
 let term1 = Complex.prod([mu1, Lambda1, besselK(0, Lambda1_b), besselK(1, Lambda2_b)]);
 let term2 = Complex.prod([mu2, Lambda2, besselK(0, Lambda2_b), besselK(1, Lambda1_b)]);
 return Complex.sub(term1, term2);
}

function scriptG(k, Lambda1, Lambda2, R, w1, mu1, mu2) {
 let b = R - w1;
 let Lambda1_b = Complex.mul(Lambda1, b);
 let Lambda2_b = Complex.mul(Lambda2, b);
 let term1 = Complex.prod([mu1, Lambda1, besselI(0, Lambda1_b), besselK(1, Lambda2_b)]);
 let term2 = Complex.prod([mu2, Lambda2, besselK(0, Lambda2_b), besselI(1, Lambda1_b)]);
 return Complex.add(term1, term2);
}

function scriptH(k, Lambda1, Lambda2, R, w1, mu1, mu2) {
 let b = R - w1;
 let Lambda1_b = Complex.mul(Lambda1, b);
 let Lambda2_b = Complex.mul(Lambda2, b);
 let term1 = Complex.prod([mu1, Lambda1, besselI(0, Lambda1_b), besselI(1, Lambda2_b)]);
 let term2 = Complex.prod([mu2, Lambda2, besselI(0, Lambda2_b), besselI(1, Lambda1_b)]);
 return Complex.sub(term1, term2);
}

function scriptI(k, Lambda1, Lambda2, R, w1, w2, mu2) {
 let c = R - (w1 + w2);
 let Lambda2_c = Complex.mul(Lambda2, c);
 let term1 = Complex.prod([mu2, Lambda2, besselK(0, Lambda2_c), besselI(1, k * c)]);
 let term2 = Complex.prod([k, besselI(0, k * c), besselK(1, Lambda2_c)]);
 return Complex.add(term1, term2);
}

function scriptK(k, Lambda1, Lambda2, R, w1, w2, mu2) {
 let c = R - (w1 + w2);
 let Lambda2_c = Complex.mul(Lambda2, c);
 let term1 = Complex.prod([mu2, Lambda2, besselI(0, Lambda2_c), besselI(1, k * c)]);
 let term2 = Complex.prod([k, besselI(0, k * c), besselI(1, Lambda2_c)]);
 return Complex.sub(term1, term2);
}


function Gamma(k, Lambda1, Lambda2, R, w1, w2, mu1, mu2) {
 let scriptE_scriptI = Complex.mul(scriptE(k, Lambda1, Lambda2, R, w1, mu1, mu2), scriptI(k, Lambda1, Lambda2, R, w1, w2, mu2));
 let scriptF_scriptK = Complex.mul(scriptF(k, Lambda1, Lambda2, R, w1, mu1, mu2), scriptK(k, Lambda1, Lambda2, R, w1, w2, mu2));
 let termEIFK = Complex.add(scriptE_scriptI, scriptF_scriptK);

 let scriptG_scriptK = Complex.mul(scriptG(k, Lambda1, Lambda2, R, w1, mu1, mu2), scriptK(k, Lambda1, Lambda2, R, w1, w2, mu2));
 let scriptH_scriptI = Complex.mul(scriptH(k, Lambda1, Lambda2, R, w1, mu1, mu2), scriptI(k, Lambda1, Lambda2, R, w1, w2, mu2));
 let termGKHI = Complex.add(scriptG_scriptK, scriptH_scriptI);

 let term1 = Complex.mul(scriptC(k, Lambda1, R, mu1), termEIFK);
 let term2 = Complex.mul(scriptD(k, Lambda1, R, mu1), termGKHI);
 let numerator = Complex.add(term1, term2)

 let term3 = Complex.mul(scriptA(k, Lambda1, R, mu1), termEIFK);
 let term4 = Complex.mul(scriptB(k, Lambda1, R, mu1), termGKHI);
 let denominator = Complex.add(term3, term4)

 return Complex.div(numerator, denominator);
}


//Axial Force
// Integrand
function axialForce(f,v,i,n,wc,l,wa,radius,w1,sigma1,mu1,w2,sigma2,mu2){
 let func = function(k) {
 let TubeFunc = Complex.sub(Gamma(k, Lambda1(k, 1, mu1, sigma1, f, v), Lambda2(k, 1, mu2, sigma2, f, v), radius, w1, w2, mu1, mu2),
 Gamma(k, Lambda1(k, -1, mu1, sigma1, f, v), Lambda2(k, -1, mu2, sigma2, f, v), radius, w1, w2, mu1, mu2)) [1];

 let CoilFunc = -mu0 * (n * i * Math.PI * Math.sin(k * l / 2) * struveL_besselK(k, radius + wa, radius + wa + wc) / k) ** 2 / (2 * k);
 return CoilFunc * TubeFunc;
 }
 // Integration
 let epsabs = 1e-5;
 let epsrel = 1e-5;
 let limit = 1000;
 let key = 1;
 let min = 0;
 // let max = Infinity;
 let max = 500; // Check that this int range makes sense.
 let intResult = [];
 let err1 = [];
 let [result, absErr] = quad(func, min, max)

 return -result;
}
window.axialForce = axialForce

//Radial Force
// Integrand
function radialForce(f, v, i, n, wc, l, wa, radius, w1, sigma1, mu1, w2, sigma2, mu2){
 let func = function(k) {
 let TubeFunc = Complex.add(Gamma(k, Lambda1(k, 1, mu1, sigma1, f, v), Lambda2(k, 1, mu2, sigma2, f, v), radius, w1, w2, mu1, mu2),
 Gamma(k, Lambda1(k, -1, mu1, sigma1, f, v), Lambda2(k, -1, mu2, sigma2, f, v), radius, w1, w2, mu1, mu2)) [0];
 let besselFunc = Complex.sub(Complex.mul((radius + wa + wc), besselK(1, k * 0.05)), Complex.mul((radius + wa), besselK(1, k * (radius + wa))));
 let CoilFunc = -mu0 * Math.PI * ((n * i * Math.sin(k * l / 2) / k) ** 2 ) * besselFunc * struveL_besselK(k, radius + wa, radius + wa + wc) / k;
 return CoilFunc * TubeFunc;
 }
 // Integration
 let epsabs = 1e-5;
 let epsrel = 1e-5;
 let limit = 1000;
 let key = 1;
 let min = 0;
 // let max = Infinity;
 let max = 500; // Check that this int range makes sense.
 let intResult = [];
 let err1 = [];
 let [result, absErr] = quad(func, min, max)

 return result;
}
window.radialForce = radialForce

// This function requires more work for the other inputs/loops.
function main(df, dv) {
 let results = [];
 for (let f of df) {
 for (let v of dv) {
 //console.log(f,v);
 results.push([R, wa, wc, mu1, sigma1, w1, mu2, sigma2, w2, f, v, Thrust(R, wa, wc, mu1, sigma1, w1, mu2, sigma2, w2, f, v)]);
 }
 }
 return results;
}
window.main = main

// let t0 = performance.now()
// // Thrust(radi, wa, wc, mu, sig, w1, mu, sig, w2, f, v)
// console.log(Thrust(0.05, 0.01, 0.01, 10, 5e6, 0.01, 25, 9e6, 0.01, 20, 5))
// let t1 = performance.now()
// console.log("Calculation took " + (t1 - t0) + " milliseconds.")

let df = [];
let dv = [];
for (let f = 0; f <= 2000; f += 200){
 df.push(f)
}
for (let v = 0; v <= 100; v += 10){
 dv.push(v)
}

// console.log(main(df,dv));

// // Note to Dan, this is another variable that I search and replace
// var $$body;
//
// global._main = main;
// module.exports = { main };

},{"../sparc-math/_es2015/sparc-math.js":95}],2:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.quad = undefined;

var _slicedToArray = function () { function sliceIterator(arr, i) { var _arr = []; var _n = true; var _d = false; var _e = undefined; try { for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) { _arr.push(_s.value); if (i && _arr.length === i) break; } } catch (err) { _d = true; _e = err; } finally { try { if (!_n && _i["return"]) _i["return"](); } finally { if (_d) throw _e; } } return _arr; } return function (arr, i) { if (Array.isArray(arr)) { return arr; } else if (Symbol.iterator in Object(arr)) { return sliceIterator(arr, i); } else { throw new TypeError("Invalid attempt to destructure non-iterable instance"); } }; }();

var _dqagse3 = require('./quadpack/dqagse.js');

var _dqagie3 = require('./quadpack/dqagie.js');

function _toConsumableArray(arr) { if (Array.isArray(arr)) { for (var i = 0, arr2 = Array(arr.length); i < arr.length; i++) { arr2[i] = arr[i]; } return arr2; } else { return Array.from(arr); } }

// Compute a definite integral.
// Integrate func from `a` to `b` (possibly infinite interval) using a
// technique from the Fortran library QUADPACK.

function quad(func, a, b) {
 var opts = arguments.length > 3 && arguments[3] !== undefined ? arguments[3] : {};

 var flip = void 0,
 defaultOpts = void 0,
 o = void 0,
 retval = void 0,
 ier = void 0,
 msg = void 0,
 msgs = void 0,
 minPoint = void 0,
 maxPoint = void 0,
 condition = void 0,
 limit = void 0;

 defaultOpts = {
 fullOutput: 0,
 epsabs: 1.49e-8,
 epsrel: 1.49e-8,
 limit: 50,
 points: null,
 weight: null,
 wvar: null,
 wopts: null,
 maxp1: 50,
 limlst: 50
 };
 o = Object.assign(defaultOpts, opts);

 var _ref = [b < a, Math.min(a, b), Math.max(a, b)];
 flip = _ref[0];
 a = _ref[1];
 b = _ref[2];


 if (o.weight === null) {
 retval = _quad(func, a, b, o.fullOutput, o.epsabs, o.epsrel, o.limit, o.points);
 } else {
 // retval = _quad_weight(func, a, b, args, fullOutput, epsabs, epsrel, limlst, limit, maxp1, weight, wvar, wopts)
 throw new Error('SciPy functionality not supported yet: weighted quad');
 }

 if (flip) retval[0] = -retval[0];
 ier = retval.pop();
 if (ier === 0) return retval;

 // Non-zero exit code, throw descriptive error
 msgs = {
 '1': 'The maximum number of subdivisions (' + o.limit + ') has been achieved.\n If increasing the limit yields no improvement it is advised to analyze\n the integrand in order to determine the difficulties. If the position of a\n local difficulty can be determined (singularity, discontinuity) one will\n probably gain from splitting up the interval and calling the integrator\n on the subranges. Perhaps a special-purpose integrator should be used.',
 '2': 'The occurrence of roundoff error is detected, which prevents\n the requested tolerance from being achieved. The error may be\n underestimated.',
 '3': 'Extremely bad integrand behavior occurs at some points of the\n integration interval.',
 '4': 'The algorithm does not converge. Roundoff error is detected\n in the extrapolation table. It is assumed that the requested tolerance\n cannot be achieved, and that the returned result (if full_output = 1) is\n the best which can be obtained.',
 '5': 'The integral is probably divergent, or slowly convergent.',
 '6': 'The input is invalid.',
 '7': 'Abnormal termination of the routine. The estimates for result\n and error are less reliable. It is assumed that the requested accuracy\n has not been achieved.',
 'unknown': 'Unknown error.'

 // TODO: implement weighted stuff
 // if weight in ['cos','sin'] and (b===Inf or a===-Inf):
 // msgs[1] = "The maximum number of cycles allowed has been achieved., e.e.\n of subintervals (a+(k-1)c, a+kc) where c = (2*int(abs(omega)+1))\n *pi/abs(omega), for k = 1, 2, ..., lst. One can allow more cycles by increasing the value of limlst. Look at info['ierlst'] with full_output=1."
 // msgs[4] = "The extrapolation table constructed for convergence acceleration\n of the series formed by the integral contributions over the cycles, \n does not converge to within the requested accuracy. Look at \n info['ierlst'] with full_output=1."
 // msgs[7] = "Bad integrand behavior occurs within one or more of the cycles.\n Location and type of the difficulty involved can be determined from \n the vector info['ierlist'] obtained with full_output=1."
 // explain = {1: "The maximum number of subdivisions (= limit) has been \n achieved on this cycle.",
 // 2: "The occurrence of roundoff error is detected and prevents\n the tolerance imposed on this cycle from being achieved.",
 // 3: "Extremely bad integrand behavior occurs at some points of\n this cycle.",
 // 4: "The integral over this cycle does not converge (to within the required accuracy) due to roundoff in the extrapolation procedure invoked on this cycle. It is assumed that the result on this interval is the best which can be obtained.",
 // 5: "The integral over this cycle is probably divergent or slowly convergent."}

 };msg = msgs[ier];
 if (msg === undefined) msg = msgs['unknown'];
 if ([1, 2, 3, 4, 5, 7].includes(ier)) {
 if (o.fullOutput) {
 // if weight in ['cos', 'sin'] and (b===Inf or a===Inf):
 // return retval[:-1] + (msg, explain)
 // else:
 return retval.concat(msg);
 } else {
 // console.warn('Integration Warning: ', msg);
 return retval;
 }
 } else if (ier === 6) {
 // # Forensic decision tree when QUADPACK throws ier=6
 if (o.epsabs <= 0) {
 // # Small error tolerance - applies to all methods
 if (o.epsrel < Math.max(50 * Number.EPSILON, 5e-29)) {
 msg = 'If \'errabs\'<=0, \'epsrel\' must be greater than both' + ('5e-29 and 50*(machine epsilon)[=' + 50 * Number.EPSILON + '].');
 }
 // elif weight in ['sin', 'cos'] and (abs(a) + abs(b)===Inf):
 // msg = ("Sine or cosine weighted intergals with infinite domain" " must have 'epsabs'>0.")
 } else if (o.weight === null) {
 if (o.points === null) {
 // QAGSE/QAGIE
 msg = 'Invalid \'limit\' argument. There must be at least one subinterval';
 } else {
 // QAGPE
 minPoint = Math.min.apply(Math, _toConsumableArray(o.points));
 maxPoint = Math.max.apply(Math, _toConsumableArray(o.points));
 condition = !(Math.min(a, b) <= minPoint && minPoint <= maxPoint && maxPoint <= Math.max(a, b));
 if (condition) {
 msg = 'All break points in \'points\' must lie within the integration limits.';
 } else if (o.points.length >= limit) {
 msg = 'Number of break points (' + o.points.length + ') must be less\n than subinterval limit (' + o.limit + ')';
 }
 }
 } else {
 if (o.maxp1 < 1) {
 msg = 'Chebyshev moment limit maxp1 must be >=1.';
 } else if (['cos', 'sin'].includes(o.weight) && Math.abs(a + b) === Infinity) {
 // QAWFE
 msg = 'Cycle limit limlst must be >=3.';
 } else if (o.weight.startswith('alg')) {
 // QAWSE
 if (Math.min(o.wvar) < -1) {
 msg = 'wvar parameters (alpha, beta) must both be >= -1.';
 }
 if (b < a) {
 msg = 'Integration limits a, b must satistfy a<b.';
 }
 } else if (o.weight === 'cauchy' && [a, b].includes(o.wvar)) {
 msg = 'Parameter \'wvar\' must not equal integration limits \'a\' or \'b\'.';
 }
 }
 }
 throw new Error(msg);
}

// RETURN: fullOutput ? [result, abserr, infoObj, ier] : [result, abserr, ier]
function _quad(func, a, b, fullOutput, epsabs, epsrel, limit, points) {
 var infbounds = void 0,
 bound = void 0,
 thePoints = void 0;

 infbounds = 0;
 if (b !== Infinity && a !== -Infinity) {
 // pass # standard integration
 } else if (b === Infinity && a !== -Infinity) {
 infbounds = 1;
 bound = a;
 } else if (b === Infinity && a === -Infinity) {
 infbounds = 2;
 bound = 0; // ignored
 } else if (b !== Infinity && a === -Infinity) {
 infbounds = -1;
 bound = b;
 } else {
 // this is the python message, pretty sure we'll never hit it but we should have
 // an error message in the else anyways
 throw new Error('Infinity comparisons don\'t work for you.');
 }

 if (points === null) {
 if (infbounds === 0) {
 return _qagse(func, a, b, fullOutput, epsabs, epsrel, limit);
 } else {
 return _qagie(func, bound, infbounds, fullOutput, epsabs, epsrel, limit);
 }
 } else {
 if (infbounds !== 0) {
 throw new Error('Infinity inputs cannot be used with break points.');
 } else {
 // Duplicates force function evaluation at singular points
 thePoints = sortUnique(points);
 thePoints = thePoints.filter(function (point) {
 return point > a;
 });
 thePoints = thePoints.filter(function (point) {
 return point < b;
 });
 thePoints = thePoints.concat([0, 0]);
 return _qagpe(func, a, b, thePoints, fullOutput, epsabs, epsrel, limit);
 }
 }
}

function _qagse(func, a, b, fullOutput, epsabs, epsrel, limit) {
 var iord = void 0,
 alist = void 0,
 blist = void 0,
 rlist = void 0,
 elist = void 0,
 infoObj = void 0;
 var neval = 0;
 var ier = 6;
 var last = 0;
 var result = 0.0;
 var abserr = 0.0;

 /* Need to check that limit is bigger than 1 */
 if (limit < 1) return [result, abserr, ier];

 /* Setup iwork and work arrays */
 iord = new Int32Array(limit);
 alist = new Float64Array(limit);
 blist = new Float64Array(limit);
 rlist = new Float64Array(limit);
 elist = new Float64Array(limit);

 var _dqagse = (0, _dqagse3.dqagse)(func, a, b, epsabs, epsrel, limit, result, abserr, neval, ier, alist, blist, rlist, elist, iord, last);

 var _dqagse2 = _slicedToArray(_dqagse, 6);

 result = _dqagse2[0];
 abserr = _dqagse2[1];
 neval = _dqagse2[2];
 ier = _dqagse2[3];
 iord = _dqagse2[4];
 last = _dqagse2[5];


 if (fullOutput) {
 infoObj = { 'neval': neval,
 'last': last,
 'iord': iord,
 'alist': alist,
 'blist': blist,
 'rlist': rlist,
 'elist': elist };
 return [result, abserr, infoObj, ier];
 } else {
 return [result, abserr, ier];
 }
}

function _qagie(func, bound, inf, fullOutput, epsabs, epsrel, limit) {
 var iord = void 0,
 alist = void 0,
 blist = void 0,
 rlist = void 0,
 elist = void 0,
 infoObj = void 0;
 var neval = 0;
 var ier = 6;
 var last = 0;
 var result = 0.0;
 var abserr = 0.0;

 /* Need to check that limit is bigger than 1 */
 if (limit < 1) return [result, abserr, ier];

 /* Setup iwork and work arrays */
 iord = new Int32Array(limit);
 alist = new Float64Array(limit);
 blist = new Float64Array(limit);
 rlist = new Float64Array(limit);
 elist = new Float64Array(limit);

 var _dqagie = (0, _dqagie3.dqagie)(func, bound, inf, epsabs, epsrel, limit, result, abserr, neval, ier, alist, blist, rlist, elist, iord, last);

 var _dqagie2 = _slicedToArray(_dqagie, 10);

 result = _dqagie2[0];
 abserr = _dqagie2[1];
 neval = _dqagie2[2];
 ier = _dqagie2[3];
 alist = _dqagie2[4];
 blist = _dqagie2[5];
 rlist = _dqagie2[6];
 elist = _dqagie2[7];
 iord = _dqagie2[8];
 last = _dqagie2[9];


 if (fullOutput) {
 infoObj = { 'neval': neval,
 'last': last,
 'iord': iord,
 'alist': alist,
 'blist': blist,
 'rlist': rlist,
 'elist': elist };
 return [result, abserr, infoObj, ier];
 } else {
 return [result, abserr, ier];
 }
}

function _qagpe(func, a, b, thePoints, fullOutput, epsabs, epsrel, limit) {
 var neval = void 0,
 ier = void 0,
 last = void 0,
 result = void 0,
 abserr = void 0,
 npts2 = void 0,
 iord = void 0,
 alist = void 0,
 blist = void 0,
 rlist = void 0,
 elist = void 0,
 pts = void 0,
 level = void 0,
 ndin = void 0,
 dqagpe = void 0;
 neval = 0;
 ier = 6;
 last = 0;
 result = 0.0;
 abserr = 0.0;

 /* Need to check that limit is bigger than 1 */
 if (limit < 1) return [result, abserr, ier];
 npts2 = thePoints.length;

 /* Setup iwork and work arrays */
 iord = new Int32Array(limit);
 alist = new Float64Array(limit);
 blist = new Float64Array(limit);
 rlist = new Float64Array(limit);
 elist = new Float64Array(limit);
 pts = new Float64Array(npts2);
 level = new Float64Array(limit);
 ndin = new Float64Array(npts2);

 var _dqagpe = dqagpe(func, a, b, npts2, thePoints, epsabs, epsrel, limit, result, abserr, neval, ier, alist, blist, rlist, elist, pts, iord, level, ndin, last);

 var _dqagpe2 = _slicedToArray(_dqagpe, 13);

 result = _dqagpe2[0];
 abserr = _dqagpe2[1];
 neval = _dqagpe2[2];
 ier = _dqagpe2[3];
 alist = _dqagpe2[4];
 blist = _dqagpe2[5];
 rlist = _dqagpe2[6];
 elist = _dqagpe2[7];
 pts = _dqagpe2[8];
 iord = _dqagpe2[9];
 level = _dqagpe2[10];
 ndin = _dqagpe2[11];
 last = _dqagpe2[12];


 if (fullOutput) {
 var infoObj = {
 'neval': neval,
 'last': last,
 'iord': iord,
 'alist': alist,
 'blist': blist,
 'rlist': rlist,
 'elist': elist,
 'pts': pts,
 'level': level,
 'ndin': ndin
 };
 return [result, abserr, infoObj, ier];
 } else {
 return [result, abserr, ier];
 }
}

// https:// stackoverflow.com/questions/4833651/javascript-array-sort-and-unique
function sortUnique(arr) {
 if (arr.length === 0) return arr;
 arr = arr.sort(function (a, b) {
 return a * 1 - b * 1;
 });
 var ret = [arr[0]];
 for (var i = 1; i < arr.length; i++) {
 // Start loop at 1: arr[0] can never be a duplicate
 if (arr[i - 1] !== arr[i]) {
 ret.push(arr[i]);
 }
 }
 return ret;
}

exports.quad = quad;
},{"./quadpack/dqagie.js":3,"./quadpack/dqagse.js":4}],3:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});

var _slicedToArray = function () { function sliceIterator(arr, i) { var _arr = []; var _n = true; var _d = false; var _e = undefined; try { for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) { _arr.push(_s.value); if (i && _arr.length === i) break; } } catch (err) { _d = true; _e = err; } finally { try { if (!_n && _i["return"]) _i["return"](); } finally { if (_d) throw _e; } } return _arr; } return function (arr, i) { if (Array.isArray(arr)) { return arr; } else if (Symbol.iterator in Object(arr)) { return sliceIterator(arr, i); } else { throw new TypeError("Invalid attempt to destructure non-iterable instance"); } }; }(); /* eslint-disable no-labels */
/* eslint-disable no-fallthrough */


exports.dqagie = dqagie;

var _d1mach = require('../../utils/fortran-utils/d1mach.js');

var _dqk15i7 = require('./dqk15i.js');

var _dqpsrt3 = require('./dqpsrt.js');

var _dqelg3 = require('./dqelg.js');

function dqagie(f, bound, inf, epsabs, epsrel, limit, result, abserr, neval, ier, alist, blist, rlist, elist, iord, last) {
 var abseps = void 0,
 area = void 0,
 area1 = void 0,
 area12 = void 0,
 area2 = void 0,
 a1 = void 0,
 a2 = void 0,
 boun = void 0,
 b1 = void 0,
 b2 = void 0,
 correc = void 0,
 defabs = void 0,
 defab1 = void 0,
 defab2 = void 0,
 dres = void 0,
 epmach = void 0,
 erlarg = void 0,
 erlast = void 0,
 errbnd = void 0,
 errmax = void 0,
 error1 = void 0,
 error2 = void 0,
 erro12 = void 0,
 errsum = void 0,
 ertest = void 0,
 oflow = void 0,
 resabs = void 0,
 reseps = void 0,
 res3la = void 0,
 rlist2 = void 0,
 small = void 0,
 uflow = void 0,
 id = void 0,
 ierro = void 0,
 iroff1 = void 0,
 iroff2 = void 0,
 iroff3 = void 0,
 jupbnd = void 0,
 k = void 0,
 ksgn = void 0,
 ktmin = void 0,
 maxerr = void 0,
 nres = void 0,
 nrmax = void 0,
 numrl2 = void 0,
 extrap = void 0,
 noext = void 0,
 condition = void 0;

 res3la = new Float64Array(3);
 rlist2 = new Float64Array(52);

 var goToLabel = 0;
 mainExecutionLoop: while (true) {
 switch (goToLabel) {
 case 0:
 epmach = (0, _d1mach.d1mach)(4);
 //
 // test on validity of parameters
 // -----------------------------
 //
 ier = 0;
 neval = 0;
 last = 0;
 result = 0.0e+00;
 abserr = 0.0e+00;
 alist[0] = 0.0e+00;
 blist[0] = 0.1e+01;
 rlist[0] = 0.0e+00;
 elist[0] = 0.0e+00;
 iord[0] = 0;
 if (epsabs <= 0.0e+00 && epsrel < Math.max(0.5e+02 * epmach, 0.5e-28)) ier = 6;
 if (ier === 6) {
 goToLabel = 999;break;
 }
 //
 //
 // first approximation to the integral
 // -----------------------------------
 //
 // determine the interval to be mapped onto (0,1).
 // if inf = 2 the integral is computed as i = i1+i2, where
 // i1 = integral of f over (-infinity,0),
 // i2 = integral of f over (0,+infinity).
 //
 boun = bound;
 if (inf === 2) boun = 0.0e+00;

 //
 // test on accuracy
 //
 var _dqk15i = (0, _dqk15i7.dqk15i)(f, boun, inf, 0.0e+00, 0.1e+01, result, abserr, defabs, resabs);

 var _dqk15i2 = _slicedToArray(_dqk15i, 4);

 result = _dqk15i2[0];
 abserr = _dqk15i2[1];
 defabs = _dqk15i2[2];
 resabs = _dqk15i2[3];
 last = 1;
 rlist[0] = result;
 elist[0] = abserr;
 iord[0] = 1;
 dres = Math.abs(result);
 errbnd = Math.max(epsabs, epsrel * dres);
 if (abserr <= 1.0e+02 * epmach * defabs && abserr > errbnd) ier = 2;
 if (limit === 1) ier = 1;
 if (ier !== 0 || abserr <= errbnd && abserr !== resabs || abserr === 0.0e+00) {
 goToLabel = 130;break;
 }
 //
 // initialization
 // --------------
 //
 uflow = (0, _d1mach.d1mach)(1);
 oflow = (0, _d1mach.d1mach)(2);
 rlist2[0] = result;
 errmax = abserr;
 maxerr = 1;
 area = result;
 errsum = abserr;
 abserr = oflow;
 nrmax = 1;
 nres = 0;
 ktmin = 0;
 numrl2 = 2;
 extrap = false;
 noext = false;
 ierro = 0;
 iroff1 = 0;
 iroff2 = 0;
 iroff3 = 0;
 ksgn = -1;
 if (dres >= (0.1e+01 - 0.5e+02 * epmach) * defabs) ksgn = 1;
 //
 // main do-loop
 // ------------
 //
 for (last = 2; last <= limit; last++) {
 //
 // bisect the subinterval with nrmax-th largest error estimate.
 //
 a1 = alist[maxerr - 1];
 b1 = 0.5e+00 * (alist[maxerr - 1] + blist[maxerr - 1]);
 a2 = b1;
 b2 = blist[maxerr - 1];
 erlast = errmax;

 var _dqk15i3 = (0, _dqk15i7.dqk15i)(f, boun, inf, a1, b1, area1, error1, resabs, defab1);

 var _dqk15i4 = _slicedToArray(_dqk15i3, 4);

 area1 = _dqk15i4[0];
 error1 = _dqk15i4[1];
 resabs = _dqk15i4[2];
 defab1 = _dqk15i4[3];

 //
 // improve previous approximations to integral
 // and error and test for accuracy.
 //
 var _dqk15i5 = (0, _dqk15i7.dqk15i)(f, boun, inf, a2, b2, area2, error2, resabs, defab2);

 var _dqk15i6 = _slicedToArray(_dqk15i5, 4);

 area2 = _dqk15i6[0];
 error2 = _dqk15i6[1];
 resabs = _dqk15i6[2];
 defab2 = _dqk15i6[3];
 area12 = area1 + area2;
 erro12 = error1 + error2;
 errsum = errsum + erro12 - errmax;
 area = area + area12 - rlist[maxerr - 1];
 if (defab1 === error1 || defab2 === error2) {
 // GO TO 15
 } else {
 condition = Math.abs(rlist[maxerr - 1] - area12) > 0.1e-04 * Math.abs(area12) || erro12 < 0.99e+00 * errmax;
 if (condition) {
 // GO TO 10
 } else {
 if (extrap) iroff2 = iroff2 + 1;
 if (!extrap) iroff1 = iroff1 + 1;
 }
 // LABEL 10
 if (last > 10 && erro12 > errmax) iroff3 = iroff3 + 1;
 }
 // LABEL 15:
 rlist[maxerr - 1] = area1;
 rlist[last - 1] = area2;
 errbnd = Math.max(epsabs, epsrel * Math.abs(area));
 //
 // test for roundoff error and eventually set error flag.
 //
 if (iroff1 + iroff2 >= 10 || iroff3 >= 20) ier = 2;
 if (iroff2 >= 5) ierro = 3;
 //
 // set error flag in the case that the number of
 // subintervals equals limit.
 //
 if (last === limit) ier = 1;
 //
 // set error flag in the case of bad integrand behaviour
 // at some points of the integration range.
 //
 condition = Math.max(Math.abs(a1), Math.abs(b2)) <= (0.1e+01 + 0.1e+03 * epmach) * (Math.abs(a2) + 0.1e+04 * uflow);
 if (condition) ier = 4;
 //
 // append the newly-created intervals to the list.
 //
 if (error2 > error1) {
 alist[maxerr - 1] = a2;
 alist[last - 1] = a1;
 blist[last - 1] = b1;
 rlist[maxerr - 1] = area2;
 rlist[last - 1] = area1;
 elist[maxerr - 1] = error2;
 elist[last - 1] = error1;
 } else {
 alist[last - 1] = a2;
 blist[maxerr - 1] = b1;
 blist[last - 1] = b2;
 elist[maxerr - 1] = error1;
 elist[last - 1] = error2;
 }
 // LABEL 30

 //
 // subroutine() dqpsrt to maintain the descending ordering
 // in the list of error estimates and select the subinterval
 // with nrmax-th largest error estimate (to be bisected next).
 //

 var _dqpsrt = (0, _dqpsrt3.dqpsrt)(limit, last, maxerr, errmax, elist, iord, nrmax);

 var _dqpsrt2 = _slicedToArray(_dqpsrt, 4);

 maxerr = _dqpsrt2[0];
 errmax = _dqpsrt2[1];
 iord = _dqpsrt2[2];
 nrmax = _dqpsrt2[3];

 if (errsum <= errbnd) {
 goToLabel = 115;break;
 }
 if (ier !== 0) {
 goToLabel = 100;break;
 }
 if (last === 2) {
 // GO TO 80
 } else {
 if (noext) continue;
 erlarg = erlarg - erlast;
 if (Math.abs(b1 - a1) > small) erlarg = erlarg + erro12;
 if (extrap) {
 // goToLabel = 40; break;
 } else {
 //
 // test whether the interval to be bisected next is the
 // smallest interval.
 //
 if (Math.abs(blist[maxerr - 1] - alist[maxerr - 1]) > small) continue;
 extrap = true;
 nrmax = 2;
 }
 if (ierro === 3 || erlarg <= ertest) {
 // goToLabel = 60; break;
 } else {
 //
 // the smallest interval has the largest error.
 // before bisecting decrease the sum of the errors over the
 // larger intervals (erlarg) and perform extrapolation.
 //
 id = nrmax;
 jupbnd = last;
 if (last > 2 + Math.trunc(limit / 2)) jupbnd = limit + 3 - last;
 for (k = id; k <= jupbnd; k++) {
 maxerr = iord[nrmax - 1];
 errmax = elist[maxerr - 1];
 if (Math.abs(blist[maxerr - 1] - alist[maxerr - 1]) > small) continue;
 nrmax = nrmax + 1;
 }
 }
 //
 // perform extrapolation.
 //
 // case 60:
 numrl2 = numrl2 + 1;
 rlist2[numrl2 - 1] = area;

 var _dqelg = (0, _dqelg3.dqelg)(numrl2, rlist2, reseps, abseps, res3la, nres);

 var _dqelg2 = _slicedToArray(_dqelg, 4);

 numrl2 = _dqelg2[0];
 reseps = _dqelg2[1];
 abseps = _dqelg2[2];
 nres = _dqelg2[3];

 ktmin = ktmin + 1;
 if (ktmin > 5 && abserr < 0.1e-02 * errsum) ier = 5;
 if (abseps >= abserr) {
 // goToLabel = 70; break;
 } else {
 ktmin = 0;
 abserr = abseps;
 result = reseps;
 correc = erlarg;
 ertest = Math.max(epsabs, epsrel * Math.abs(reseps));
 if (abserr <= ertest) {
 goToLabel = 100;break;
 }
 }
 //
 // prepare bisection of the smallest interval.
 //
 if (numrl2 === 1) noext = true;
 if (ier === 5) {
 goToLabel = 100;break;
 }
 maxerr = iord[0];
 errmax = elist[maxerr - 1];
 nrmax = 1;
 extrap = false;
 small = small * 0.5e+00;
 erlarg = errsum;
 continue;
 }
 // case 80:
 small = 0.375e+00;
 erlarg = errsum;
 ertest = errbnd;
 rlist2[1] = area;
 } // end of loop
 if (goToLabel > 100) break;
 case 100:
 //
 // set final result and error estimate.
 // ------------------------------------
 //
 if (abserr === oflow) {
 goToLabel = 115;break;
 }
 if (ier + ierro === 0) {
 goToLabel = 110;break;
 }
 if (ierro === 3) abserr = abserr + correc;
 if (ier === 0) ier = 3;
 if (result !== 0.0e+00 && area !== 0.0e+00) {
 goToLabel = 105;break;
 }
 if (abserr > errsum) {
 goToLabel = 115;break;
 }
 if (area === 0.0e+00) {
 goToLabel = 130;break;
 }
 goToLabel = 110;break;
 case 105:
 if (abserr / Math.abs(result) > errsum / Math.abs(area)) {
 goToLabel = 115;break;
 }
 //
 // test on divergence
 //
 case 110:
 if (ksgn === -1 && Math.max(Math.abs(result), Math.abs(area)) <= defabs * 0.1e-01) {
 goToLabel = 130;break;
 }
 if (result / area < 0.1e-01 || result / area > 0.1e+03 || errsum > Math.abs(area)) ier = 6;
 goToLabel = 130;break;
 //
 // compute global integral sum.
 //
 case 115:
 result = 0.0e+00;
 for (k = 1; k <= last; k++) {
 result = result + rlist[k - 1];
 }
 abserr = errsum;
 case 130:
 neval = 30 * last - 15;
 if (inf === 2) neval = 2 * neval;
 if (ier > 2) ier = ier - 1;
 case 999:

 default:
 break mainExecutionLoop;
 }
 }
 return [result, abserr, neval, ier, alist, blist, rlist, elist, iord, last];
}
},{"../../utils/fortran-utils/d1mach.js":91,"./dqelg.js":5,"./dqk15i.js":6,"./dqpsrt.js":8}],4:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});

var _slicedToArray = function () { function sliceIterator(arr, i) { var _arr = []; var _n = true; var _d = false; var _e = undefined; try { for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) { _arr.push(_s.value); if (i && _arr.length === i) break; } } catch (err) { _d = true; _e = err; } finally { try { if (!_n && _i["return"]) _i["return"](); } finally { if (_d) throw _e; } } return _arr; } return function (arr, i) { if (Array.isArray(arr)) { return arr; } else if (Symbol.iterator in Object(arr)) { return sliceIterator(arr, i); } else { throw new TypeError("Invalid attempt to destructure non-iterable instance"); } }; }(); /* eslint-disable no-labels */
/* eslint-disable no-fallthrough */


exports.dqagse = dqagse;

var _d1mach = require('../../utils/fortran-utils/d1mach.js');

var _dqk7 = require('./dqk21.js');

var _dqpsrt3 = require('./dqpsrt.js');

var _dqelg3 = require('./dqelg.js');

function dqagse(f, a, b, epsabs, epsrel, limit, result, abserr, neval, ier, alist, blist, rlist, elist, iord, last) {
 var abseps = void 0,
 area = void 0,
 area1 = void 0,
 area12 = void 0,
 area2 = void 0,
 a1 = void 0,
 a2 = void 0,
 b1 = void 0,
 b2 = void 0,
 correc = void 0,
 defabs = void 0,
 defab1 = void 0,
 defab2 = void 0,
 dres = void 0,
 epmach = void 0,
 erlarg = void 0,
 erlast = void 0,
 errbnd = void 0,
 errmax = void 0,
 error1 = void 0,
 error2 = void 0,
 erro12 = void 0,
 errsum = void 0,
 ertest = void 0,
 oflow = void 0,
 resabs = void 0,
 reseps = void 0,
 res3la = void 0,
 rlist2 = void 0,
 small = void 0,
 uflow = void 0,
 id = void 0,
 ierro = void 0,
 iroff1 = void 0,
 iroff2 = void 0,
 iroff3 = void 0,
 jupbnd = void 0,
 k = void 0,
 ksgn = void 0,
 ktmin = void 0,
 maxerr = void 0,
 nres = void 0,
 nrmax = void 0,
 numrl2 = void 0,
 extrap = void 0,
 noext = void 0;

 res3la = new Float64Array(3);
 rlist2 = new Float64Array(52);

 var goToLabel = 0;
 mainExecutionLoop: while (true) {
 switch (goToLabel) {
 case 0:
 epmach = (0, _d1mach.d1mach)(4);
 //
 // test on validity of parameters
 // ------------------------------
 ier = 0;
 neval = 0;
 last = 0;
 result = 0.0e+00;
 abserr = 0.0e+00;
 alist[0] = a;
 blist[0] = b;
 rlist[0] = 0.0e+00;
 elist[0] = 0.0e+00;
 if (epsabs <= 0.0e+00 && epsrel < Math.max(0.5e+02 * epmach, 0.5e-28)) ier = 6;
 if (ier === 6) {
 goToLabel = 999;break;
 }
 //
 // first approximation to the integral
 // -----------------------------------
 //
 uflow = (0, _d1mach.d1mach)(1);
 oflow = (0, _d1mach.d1mach)(2);
 ierro = 0;

 //
 // test on accuracy.
 //
 var _dqk = (0, _dqk7.dqk21)(f, a, b, result, abserr, defabs, resabs);

 var _dqk2 = _slicedToArray(_dqk, 4);

 result = _dqk2[0];
 abserr = _dqk2[1];
 defabs = _dqk2[2];
 resabs = _dqk2[3];
 dres = Math.abs(result);
 errbnd = Math.max(epsabs, epsrel * dres);
 last = 1;
 rlist[0] = result;
 elist[0] = abserr;
 iord[0] = 1;
 if (abserr <= 1.0e+02 * epmach * defabs && abserr > errbnd) ier = 2;
 if (limit === 1) ier = 1;
 if (ier !== 0 || abserr <= errbnd && abserr !== resabs || abserr === 0) {
 goToLabel = 140;break;
 }
 //
 // initialization
 // --------------
 //
 rlist2[0] = result;
 errmax = abserr;
 maxerr = 1;
 area = result;
 errsum = abserr;
 abserr = oflow;
 nrmax = 1;
 nres = 0;
 numrl2 = 2;
 ktmin = 0;
 extrap = false;
 noext = false;
 iroff1 = 0;
 iroff2 = 0;
 iroff3 = 0;
 ksgn = -1;
 if (dres >= (0.1e+01 - 0.5e+02 * epmach) * defabs) ksgn = 1;
 //
 // main do-loop
 // ------------
 //
 mainDoLoop: for (last = 2; last <= limit; last++) {
 //
 // bisect the subinterval with the nrmax-th largest error
 // estimate.
 //
 a1 = alist[maxerr - 1];
 b1 = 0.5e+00 * (alist[maxerr - 1] + blist[maxerr - 1]);
 a2 = b1;
 b2 = blist[maxerr - 1];
 erlast = errmax;

 var _dqk3 = (0, _dqk7.dqk21)(f, a1, b1, area1, error1, resabs, defab1);

 var _dqk4 = _slicedToArray(_dqk3, 4);

 area1 = _dqk4[0];
 error1 = _dqk4[1];
 resabs = _dqk4[2];
 defab1 = _dqk4[3];

 //
 // improve previous approximations to integral
 // and error and test for accuracy.
 //
 var _dqk5 = (0, _dqk7.dqk21)(f, a2, b2, area2, error2, resabs, defab2);

 var _dqk6 = _slicedToArray(_dqk5, 4);

 area2 = _dqk6[0];
 error2 = _dqk6[1];
 resabs = _dqk6[2];
 defab2 = _dqk6[3];
 area12 = area1 + area2;
 erro12 = error1 + error2;
 errsum = errsum + erro12 - errmax;
 area = area + area12 - rlist[maxerr - 1];
 if (defab1 === error1 || defab2 === error2) {
 // goToLabel = 15; break;
 } else {
 if (Math.abs(rlist[maxerr - 1] - area12) > 0.1e-04 * Math.abs(area12) || erro12 < 0.99e+00 * errmax) {
 // goToLabel = 10; break;
 } else {
 if (extrap) iroff2 = iroff2 + 1;
 if (!extrap) iroff1 = iroff1 + 1;
 }
 // case 10:
 if (last > 10 && erro12 > errmax) iroff3 = iroff3 + 1;
 }
 // case 15:
 rlist[maxerr - 1] = area1;
 rlist[last - 1] = area2;
 errbnd = Math.max(epsabs, epsrel * Math.abs(area));
 //
 // test for roundoff error and eventually set error flag.
 //
 if (iroff1 + iroff2 >= 10 || iroff3 >= 20) ier = 2;
 if (iroff2 >= 5) ierro = 3;
 //
 // set error flag in the case that the number of subintervals
 // equals limit.
 //
 if (last === limit) ier = 1;
 //
 // set error flag in the case of bad integrand behaviour
 // at a point of the integration range.
 //
 if (Math.max(Math.abs(a1), Math.abs(b2)) <= (0.1e+01 + 0.1e+03 * epmach) * (Math.abs(a2) + 0.1e+04 * uflow)) ier = 4;
 //
 // append the newly-created intervals to the list.
 //
 if (error2 > error1) {
 alist[maxerr - 1] = a2;
 alist[last - 1] = a1;
 blist[last - 1] = b1;
 rlist[maxerr - 1] = area2;
 rlist[last - 1] = area1;
 elist[maxerr - 1] = error2;
 elist[last - 1] = error1;
 } else {
 alist[last - 1] = a2;
 blist[maxerr - 1] = b1;
 blist[last - 1] = b2;
 elist[maxerr - 1] = error1;
 elist[last - 1] = error2;
 }
 //
 // subroutine() dqpsrt to maintain the descending ordering
 // in the list of error estimates and select the subinterval
 // with nrmax-th largest error estimate (to be bisected next).
 //

 // ***jump out of do-loop
 var _dqpsrt = (0, _dqpsrt3.dqpsrt)(limit, last, maxerr, errmax, elist, iord, nrmax);

 var _dqpsrt2 = _slicedToArray(_dqpsrt, 4);

 maxerr = _dqpsrt2[0];
 errmax = _dqpsrt2[1];
 iord = _dqpsrt2[2];
 nrmax = _dqpsrt2[3];
 if (errsum <= errbnd) {
 goToLabel = 115;break;
 }
 // ***jump out of do-loop
 if (ier !== 0) {
 goToLabel = 100;break;
 }
 if (last === 2) {
 // goToLabel = 80; break;
 } else {
 if (noext) continue;
 erlarg = erlarg - erlast;
 if (Math.abs(b1 - a1) > small) erlarg = erlarg + erro12;
 if (extrap) {
 // goToLabel = 40; break;
 } else {
 //
 // test whether the interval to be bisected next is the
 // smallest interval.
 //
 if (Math.abs(blist[maxerr - 1] - alist[maxerr - 1]) > small) continue;
 extrap = true;
 nrmax = 2;
 }
 // case 40:
 if (ierro === 3 || erlarg <= ertest) {
 // goToLabel = 60; break;
 } else {
 //
 // the smallest interval has the largest error.
 // before bisecting decrease the sum of the errors over the
 // larger intervals (erlarg) and perform extrapolation.
 //
 id = nrmax;
 jupbnd = last;
 if (last > 2 + Math.trunc(limit / 2)) jupbnd = limit + 3 - last;
 for (k = id; k <= jupbnd; k++) {
 maxerr = iord[nrmax - 1];
 errmax = elist[maxerr - 1];
 // ***jump out of do-loop
 if (Math.abs(blist[maxerr - 1] - alist[maxerr - 1]) > small) continue mainDoLoop;
 nrmax = nrmax + 1;
 }
 }

 //
 // perform extrapolation.
 //
 // case 60:
 numrl2 = numrl2 + 1;
 rlist2[numrl2 - 1] = area;

 var _dqelg = (0, _dqelg3.dqelg)(numrl2, rlist2, reseps, abseps, res3la, nres);

 var _dqelg2 = _slicedToArray(_dqelg, 4);

 numrl2 = _dqelg2[0];
 reseps = _dqelg2[1];
 abseps = _dqelg2[2];
 nres = _dqelg2[3];

 ktmin = ktmin + 1;
 if (ktmin > 5 && abserr < 0.1e-02 * errsum) ier = 5;
 if (abseps >= abserr) {
 // goToLabel = 70; break;
 } else {
 ktmin = 0;
 abserr = abseps;
 result = reseps;
 correc = erlarg;
 ertest = Math.max(epsabs, epsrel * Math.abs(reseps));
 // ***jump out of do-loop
 if (abserr <= ertest) {
 goToLabel = 100;break;
 }
 }

 //
 // prepare bisection of the smallest interval.
 //
 // case 70:
 if (numrl2 === 1) noext = true;
 if (ier === 5) {
 goToLabel = 100;break;
 }
 maxerr = iord[0];
 errmax = elist[maxerr - 1];
 nrmax = 1;
 extrap = false;
 small = small * 0.5e+00;
 erlarg = errsum;
 continue;
 }
 // case 80:
 small = Math.abs(b - a) * 0.375e+00;
 erlarg = errsum;
 ertest = errbnd;
 rlist2[1] = area;
 }
 if (goToLabel > 100) break;
 //
 // set final result and error estimate.
 // ------------------------------------
 //
 case 100:
 if (abserr === oflow) {
 goToLabel = 115;break;
 }
 if (ier + ierro === 0) {
 goToLabel = 110;break;
 }
 if (ierro === 3) abserr = abserr + correc;
 if (ier === 0) ier = 3;
 if (result !== 0.0e+00 && area !== 0.0e+00) {
 goToLabel = 105;break;
 }
 if (abserr > errsum) {
 goToLabel = 115;break;
 }
 if (area === 0.0e+00) {
 goToLabel = 130;break;
 }
 goToLabel = 110;break;
 case 105:
 if (abserr / Math.abs(result) > errsum / Math.abs(area)) {
 goToLabel = 115;break;
 }
 //
 // test on divergence.
 //
 case 110:
 if (ksgn === -1 && Math.max(Math.abs(result), Math.abs(area)) <= defabs * 0.1e-01) {
 goToLabel = 130;break;
 }
 if (result / area < 0.1e-01 || result / area > 0.1e+03 || errsum > Math.abs(area)) ier = 6;
 goToLabel = 130;break;
 //
 // compute global integral sum.
 //
 case 115:
 result = 0.0e+00;
 for (k = 1; k <= last; k++) {
 result = result + rlist[k - 1];
 }
 abserr = errsum;
 case 130:
 if (ier > 2) ier = ier - 1;
 case 140:
 neval = 42 * last - 21;
 case 999:
 default:
 break mainExecutionLoop;
 }
 }
 return [result, abserr, neval, ier, iord, last];
}
},{"../../utils/fortran-utils/d1mach.js":91,"./dqelg.js":5,"./dqk21.js":7,"./dqpsrt.js":8}],5:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.dqelg = dqelg;

var _d1mach = require('../../utils/fortran-utils/d1mach.js');

function dqelg(n, epstab, result, abserr, res3la, nres) {
 var delta1 = void 0,
 delta2 = void 0,
 delta3 = void 0,
 epmach = void 0,
 epsinf = void 0,
 error = void 0,
 err1 = void 0,
 err2 = void 0,
 err3 = void 0,
 e0 = void 0,
 e1 = void 0,
 e1abs = void 0,
 e2 = void 0,
 e3 = void 0,
 oflow = void 0,
 res = void 0,
 ss = void 0,
 tol1 = void 0,
 tol2 = void 0,
 tol3 = void 0,
 i = void 0,
 ib = void 0,
 ib2 = void 0,
 ie = void 0,
 indx = void 0,
 k1 = void 0,
 k2 = void 0,
 k3 = void 0,
 limexp = void 0,
 newelm = void 0,
 num = void 0;
 var goToLabel = 0;
 mainExecutionLoop: while (true) {
 switch (goToLabel) {
 case 0:
 epmach = (0, _d1mach.d1mach)(4);
 oflow = (0, _d1mach.d1mach)(2);
 nres = nres + 1;
 abserr = oflow;
 result = epstab[n - 1];
 if (n < 3) {
 goToLabel = 100;break;
 }
 limexp = 50;
 epstab[n + 1] = epstab[n - 1];
 newelm = (n - 1) / 2;
 epstab[n - 1] = oflow;
 num = n;
 k1 = n;
 for (i = 1; i <= newelm; i++) {
 k2 = k1 - 1;
 k3 = k1 - 2;
 res = epstab[k1 + 1];
 e0 = epstab[k3 - 1];
 e1 = epstab[k2 - 1];
 e2 = res;
 e1abs = Math.abs(e1);
 delta2 = e2 - e1;
 err2 = Math.abs(delta2);
 tol2 = Math.max(Math.abs(e2), e1abs) * epmach;
 delta3 = e1 - e0;
 err3 = Math.abs(delta3);
 tol3 = Math.max(e1abs, Math.abs(e0)) * epmach;
 if (err2 > tol2 || err3 > tol3) {
 // goToLabel = 10; break;
 } else {
 //
 // if e0, e1 and e2 are equal to within machine
 // accuracy, convergence is assumed.
 // result = e2
 // abserr = Math.abs(e1-e0)+Math.abs(e2-e1)
 //
 result = res;
 abserr = err2 + err3;
 // ***jump out of do-loop
 goToLabel = 100;break;
 }
 // case 10:
 e3 = epstab[k1 - 1];
 epstab[k1 - 1] = e1;
 delta1 = e1 - e3;
 err1 = Math.abs(delta1);
 tol1 = Math.max(e1abs, Math.abs(e3)) * epmach;
 //
 // if two elements are very close to each other, omit
 // a part of the table by adjusting the value of n
 //
 if (err1 <= tol1 || err2 <= tol2 || err3 <= tol3) {
 n = i + i - 1;
 // ***jump out of do-loop
 goToLabel = 50;break;
 }
 ss = 0.1e+01 / delta1 + 0.1e+01 / delta2 - 0.1e+01 / delta3;
 epsinf = Math.abs(ss * e1);
 //
 // test to detect irregular behaviour in the table, and
 // eventually omit a part of the table adjusting the value
 // of n.
 //
 if (epsinf > 0.1e-03) {
 // goToLabel = 30; break;
 } else {
 n = i + i - 1;
 // ***jump out of do-loop
 goToLabel = 50;break;
 }
 //
 // compute a new element and eventually adjust
 // the value of result.
 //
 // case 30:
 res = e1 + 0.1e+01 / ss;
 epstab[k1 - 1] = res;
 k1 = k1 - 2;
 error = err2 + Math.abs(res - e2) + err3;
 if (error > abserr) continue;
 abserr = error;
 result = res;
 }
 if (goToLabel > 50) break;
 //
 // shift the table.
 //
 case 50:
 if (n === limexp) n = 2 * Math.trunc(limexp / 2) - 1;
 ib = 1;
 if (Math.trunc(num / 2) * 2 === num) ib = 2;
 ie = newelm + 1;
 for (i = 1; i <= ie; i++) {
 ib2 = ib + 2;
 epstab[ib - 1] = epstab[ib2 - 1];
 ib = ib2;
 }
 if (num === n) {
 goToLabel = 80;break;
 }
 indx = num - n + 1;
 for (i = 1; i <= n; i++) {
 epstab[i - 1] = epstab[indx - 1];
 indx = indx + 1;
 }
 case 80:
 if (nres >= 4) {
 goToLabel = 90;break;
 }
 res3la[nres - 1] = result;
 abserr = oflow;
 goToLabel = 100;break;
 //
 // compute error estimate
 //
 case 90:
 abserr = Math.abs(result - res3la[2]) + Math.abs(result - res3la[1]) + Math.abs(result - res3la[0]);
 res3la[0] = res3la[1];
 res3la[1] = res3la[2];
 res3la[2] = result;
 case 100:
 abserr = Math.max(abserr, 0.5e+01 * epmach * Math.abs(result));
 default:
 break mainExecutionLoop;
 }
 }
 return [n, result, abserr, nres];
} /* eslint-disable no-labels */
/* eslint-disable no-fallthrough */
},{"../../utils/fortran-utils/d1mach.js":91}],6:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.dqk15i = dqk15i;

var _d1mach = require('../../utils/fortran-utils/d1mach.js');

function dqk15i(f, boun, inf, a, b, result, abserr, resabs, resasc) {
 var absc = void 0,
 absc1 = void 0,
 absc2 = void 0,
 centr = void 0,
 dinf = void 0,
 epmach = void 0,
 fc = void 0,
 fsum = void 0,
 fval1 = void 0,
 fval2 = void 0,
 fv1 = void 0,
 fv2 = void 0,
 hlgth = void 0,
 resg = void 0,
 resk = void 0,
 reskh = void 0,
 tabsc1 = void 0,
 tabsc2 = void 0,
 uflow = void 0,
 wg = void 0,
 wgk = void 0,
 xgk = void 0,
 j = void 0;

 fv1 = new Float64Array(7);
 fv2 = new Float64Array(7);
 xgk = new Float64Array(8);
 wgk = new Float64Array(8);
 wg = new Float64Array(8);

 var goToLabel = 0;
 mainExecutionLoop: while (true) {
 switch (goToLabel) {
 case 0:
 wg[0] = 0.0e0;
 wg[1] = 0.129484966168869693270611432679082e0;
 wg[2] = 0.0e0;
 wg[3] = 0.279705391489276667901467771423780e0;
 wg[4] = 0.0e0;
 wg[5] = 0.381830050505118944950369775488975e0;
 wg[6] = 0.0e0;
 wg[7] = 0.417959183673469387755102040816327e0;
 //
 xgk[0] = 0.991455371120812639206854697526329e0;
 xgk[1] = 0.949107912342758524526189684047851e0;
 xgk[2] = 0.864864423359769072789712788640926e0;
 xgk[3] = 0.741531185599394439863864773280788e0;
 xgk[4] = 0.586087235467691130294144838258730e0;
 xgk[5] = 0.405845151377397166906606412076961e0;
 xgk[6] = 0.207784955007898467600689403773245e0;
 xgk[7] = 0.000000000000000000000000000000000e0;
 //
 wgk[0] = 0.022935322010529224963732008058970e0;
 wgk[1] = 0.063092092629978553290700663189204e0;
 wgk[2] = 0.104790010322250183839876322541518e0;
 wgk[3] = 0.140653259715525918745189590510238e0;
 wgk[4] = 0.169004726639267902826583426598550e0;
 wgk[5] = 0.190350578064785409913256402421014e0;
 wgk[6] = 0.204432940075298892414161999234649e0;
 wgk[7] = 0.209482141084727828012999174891714e0;
 //
 epmach = (0, _d1mach.d1mach)(4);
 uflow = (0, _d1mach.d1mach)(1);
 dinf = Math.min(1, inf);
 //
 centr = 0.5e+00 * (a + b);
 hlgth = 0.5e+00 * (b - a);
 tabsc1 = boun + dinf * (0.1e+01 - centr) / centr;
 fval1 = f(tabsc1);
 if (inf === 2) fval1 = fval1 + f(-tabsc1);
 fc = fval1 / centr / centr;
 //
 // compute the 15-point kronrod approximation to
 // the integral, and estimate the error.
 //
 resg = wg[7] * fc;
 resk = wgk[7] * fc;
 resabs = Math.abs(resk);
 for (j = 1; j <= 7; j++) {
 absc = hlgth * xgk[j - 1];
 absc1 = centr - absc;
 absc2 = centr + absc;
 tabsc1 = boun + dinf * (0.1e+01 - absc1) / absc1;
 tabsc2 = boun + dinf * (0.1e+01 - absc2) / absc2;
 fval1 = f(tabsc1);
 fval2 = f(tabsc2);
 if (isNaN(fval1)) {}
 if (inf === 2) fval1 = fval1 + f(-tabsc1);
 if (inf === 2) fval2 = fval2 + f(-tabsc2);
 fval1 = fval1 / absc1 / absc1;
 fval2 = fval2 / absc2 / absc2;
 fv1[j - 1] = fval1;
 fv2[j - 1] = fval2;
 fsum = fval1 + fval2;
 resg = resg + wg[j - 1] * fsum;
 resk = resk + wgk[j - 1] * fsum;
 resabs = resabs + wgk[j - 1] * (Math.abs(fval1) + Math.abs(fval2));
 }
 reskh = resk * 0.5e+00;
 resasc = wgk[7] * Math.abs(fc - reskh);
 for (j = 1; j <= 7; j++) {
 resasc = resasc + wgk[j - 1] * (Math.abs(fv1[j - 1] - reskh) + Math.abs(fv2[j - 1] - reskh));
 }
 result = resk * hlgth;
 resasc = resasc * hlgth;
 resabs = resabs * hlgth;
 abserr = Math.abs((resk - resg) * hlgth);
 if (resasc !== 0.0e+00 && abserr !== 0.e0) {
 abserr = resasc * Math.min(0.1e+01, (0.2e+03 * abserr / resasc) ** 1.5e+00);
 }
 if (resabs > uflow / (0.5e+02 * epmach)) abserr = Math.max(epmach * 0.5e+02 * resabs, abserr);
 default:
 break mainExecutionLoop;
 }
 }
 return [result, abserr, resabs, resasc];
} /* eslint-disable no-labels */
/* eslint-disable no-fallthrough */
},{"../../utils/fortran-utils/d1mach.js":91}],7:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.dqk21 = dqk21;

var _d1mach = require('../../utils/fortran-utils/d1mach.js');

function dqk21(f, a, b, result, abserr, resabs, resasc) {
 var absc = void 0,
 centr = void 0,
 dhlgth = void 0,
 epmach = void 0,
 fc = void 0,
 fsum = void 0,
 fval1 = void 0,
 fval2 = void 0,
 fv1 = void 0,
 fv2 = void 0,
 hlgth = void 0,
 resg = void 0,
 resk = void 0,
 reskh = void 0,
 uflow = void 0,
 wg = void 0,
 wgk = void 0,
 xgk = void 0,
 j = void 0,
 jtw = void 0,
 jtwm1 = void 0;

 fv1 = new Float64Array(10);
 fv2 = new Float64Array(10);
 wg = new Float64Array(5);
 wgk = new Float64Array(11);
 xgk = new Float64Array(11);

 var goToLabel = 0;
 mainExecutionLoop: while (true) {
 switch (goToLabel) {
 case 0:
 wg[0] = 0.066671344308688137593568809893332e0;
 wg[1] = 0.149451349150580593145776339657697e0;
 wg[2] = 0.219086362515982043995534934228163e0;
 wg[3] = 0.269266719309996355091226921569469e0;
 wg[4] = 0.295524224714752870173892994651338e0;
 //
 xgk[0] = 0.995657163025808080735527280689003e0;
 xgk[1] = 0.973906528517171720077964012084452e0;
 xgk[2] = 0.930157491355708226001207180059508e0;
 xgk[3] = 0.865063366688984510732096688423493e0;
 xgk[4] = 0.780817726586416897063717578345042e0;
 xgk[5] = 0.679409568299024406234327365114874e0;
 xgk[6] = 0.562757134668604683339000099272694e0;
 xgk[7] = 0.433395394129247190799265943165784e0;
 xgk[8] = 0.294392862701460198131126603103866e0;
 xgk[9] = 0.148874338981631210884826001129720e0;
 xgk[10] = 0.000000000000000000000000000000000e0;
 //
 wgk[0] = 0.011694638867371874278064396062192e0;
 wgk[1] = 0.032558162307964727478818972459390e0;
 wgk[2] = 0.054755896574351996031381300244580e0;
 wgk[3] = 0.075039674810919952767043140916190e0;
 wgk[4] = 0.093125454583697605535065465083366e0;
 wgk[5] = 0.109387158802297641899210590325805e0;
 wgk[6] = 0.123491976262065851077958109831074e0;
 wgk[7] = 0.134709217311473325928054001771707e0;
 wgk[8] = 0.142775938577060080797094273138717e0;
 wgk[9] = 0.147739104901338491374841515972068e0;
 wgk[10] = 0.149445554002916905664936468389821e0;
 //
 epmach = (0, _d1mach.d1mach)(4);
 uflow = (0, _d1mach.d1mach)(1);
 //
 centr = 0.5e+00 * (a + b);
 hlgth = 0.5e+00 * (b - a);
 dhlgth = Math.abs(hlgth);
 //
 // compute the 21-point kronrod approximation to
 // the integral, and estimate the absolute error.
 //
 resg = 0.0e+00;
 fc = f(centr);
 resk = wgk[10] * fc;
 resabs = Math.abs(resk);
 for (j = 1; j <= 5; j++) {
 jtw = 2 * j;
 absc = hlgth * xgk[jtw - 1];
 fval1 = f(centr - absc);
 fval2 = f(centr + absc);
 fv1[jtw - 1] = fval1;
 fv2[jtw - 1] = fval2;
 fsum = fval1 + fval2;
 resg = resg + wg[j - 1] * fsum;
 resk = resk + wgk[jtw - 1] * fsum;
 resabs = resabs + wgk[jtw - 1] * (Math.abs(fval1) + Math.abs(fval2));
 }
 for (j = 1; j <= 5; j++) {
 jtwm1 = 2 * j - 1;
 absc = hlgth * xgk[jtwm1 - 1];
 fval1 = f(centr - absc);
 fval2 = f(centr + absc);
 fv1[jtwm1 - 1] = fval1;
 fv2[jtwm1 - 1] = fval2;
 fsum = fval1 + fval2;
 resk = resk + wgk[jtwm1 - 1] * fsum;
 resabs = resabs + wgk[jtwm1 - 1] * (Math.abs(fval1) + Math.abs(fval2));
 }
 reskh = resk * 0.5e+00;
 resasc = wgk[10] * Math.abs(fc - reskh);
 for (j = 1; j <= 10; j++) {
 resasc = resasc + wgk[j - 1] * (Math.abs(fv1[j - 1] - reskh) + Math.abs(fv2[j - 1] - reskh));
 }
 result = resk * hlgth;
 resabs = resabs * dhlgth;
 resasc = resasc * dhlgth;
 abserr = Math.abs((resk - resg) * hlgth);
 if (resasc !== 0.0e+00 && abserr !== 0.0e+00) {
 abserr = resasc * Math.min(0.1e+01, (0.2e+03 * abserr / resasc) ** 1.5e+00);
 }
 if (resabs > uflow / (0.5e+02 * epmach)) abserr = Math.max(epmach * 0.5e+02 * resabs, abserr);
 default:
 break mainExecutionLoop;
 }
 }
 return [result, abserr, resabs, resasc];
} /* eslint-disable no-labels */
/* eslint-disable no-fallthrough */
},{"../../utils/fortran-utils/d1mach.js":91}],8:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.dqpsrt = dqpsrt;
/* eslint-disable no-labels */
/* eslint-disable no-fallthrough */
function dqpsrt(limit, last, maxerr, ermax, elist, iord, nrmax) {
 var errmax = void 0,
 errmin = void 0,
 i = void 0,
 ibeg = void 0,
 ido = void 0,
 isucc = void 0,
 j = void 0,
 jbnd = void 0,
 jupbn = void 0,
 k = void 0;
 var goToLabel = 0;
 mainExecutionLoop: while (true) {
 switch (goToLabel) {
 case 0:
 //
 // check whether the list contains more than
 // two error estimates.
 //
 if (last > 2) {
 goToLabel = 10;break;
 }
 iord[0] = 1;
 iord[1] = 2;
 goToLabel = 90;break;
 //
 // this part of the routine is only executed if, due to a
 // difficult integrand, subdivision increased the error
 // estimate. in the normal case the insert procedure should
 // start after the nrmax-th largest error estimate.
 //
 case 10:
 errmax = elist[maxerr - 1];
 if (nrmax === 1) {
 goToLabel = 30;break;
 }
 ido = nrmax - 1;
 for (i = 1; i <= ido; i++) {
 isucc = iord[nrmax - 2];
 // ***jump out of do-loop
 if (errmax <= elist[isucc - 1]) {
 goToLabel = 30;break;
 }
 iord[nrmax - 1] = isucc;
 nrmax = nrmax - 1;
 }

 //
 // compute the number of elements in the list to be maintained
 // in descending order. this number depends on the number of
 // subdivisions still allowed.
 //
 case 30:
 jupbn = last;
 if (last > Math.trunc(limit / 2) + 2) jupbn = limit + 3 - last;
 errmin = elist[last - 1];
 //
 // insert errmax by traversing the list top-down,
 // starting comparison from the element elist(iord(nrmax+1)).
 //
 jbnd = jupbn - 1;
 ibeg = nrmax + 1;
 if (ibeg > jbnd) {
 goToLabel = 50;break;
 }
 for (i = ibeg; i <= jbnd; i++) {
 isucc = iord[i - 1];
 // ***jump out of do-loop
 if (errmax >= elist[isucc - 1]) {
 goToLabel = 60;break;
 }
 iord[i - 2] = isucc;
 }
 if (goToLabel === 60) break;
 case 50:
 iord[jbnd - 1] = maxerr;
 iord[jupbn - 1] = last;
 goToLabel = 90;break;
 //
 // insert errmin by traversing the list bottom-up.
 //
 case 60:
 iord[i - 2] = maxerr;
 k = jbnd;
 for (j = i; j <= jbnd; j++) {
 isucc = iord[k - 1];
 // ***jump out of do-loop
 if (errmin < elist[isucc - 1]) {
 goToLabel = 80;break;
 }
 iord[k] = isucc;
 k = k - 1;
 }
 if (goToLabel === 80) break;
 iord[i - 1] = last;
 goToLabel = 90;break;
 case 80:
 iord[k] = last;
 //
 // set maxerr and ermax.
 //
 case 90:
 maxerr = iord[nrmax - 1];
 ermax = elist[maxerr - 1];
 default:
 break mainExecutionLoop;
 }
 }
 return [maxerr, ermax, iord, nrmax];
}
},{}],9:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.airyB = exports.airyA = exports.besselH = exports.besselK = exports.besselI = exports.besselY = exports.besselJ = undefined;

var _slicedToArray = function () { function sliceIterator(arr, i) { var _arr = []; var _n = true; var _d = false; var _e = undefined; try { for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) { _arr.push(_s.value); if (i && _arr.length === i) break; } } catch (err) { _d = true; _e = err; } finally { try { if (!_n && _i["return"]) _i["return"](); } finally { if (_d) throw _e; } } return _arr; } return function (arr, i) { if (Array.isArray(arr)) { return arr; } else if (Symbol.iterator in Object(arr)) { return sliceIterator(arr, i); } else { throw new TypeError("Invalid attempt to destructure non-iterable instance"); } }; }();

var _complex = require('../utils/complex.js');

var _complex2 = _interopRequireDefault(_complex);

var _zbesj = require('./amos/zbesj.js');

var _zbesy = require('./amos/zbesy.js');

var _zbesi = require('./amos/zbesi.js');

var _zbesk = require('./amos/zbesk.js');

var _zbesh = require('./amos/zbesh.js');

var _zairy = require('./amos/zairy.js');

var _zbiry = require('./amos/zbiry.js');

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

var INPUT_ERROR = 'INPUT ERROR - NO COMPUTATION';
var OVERFLOW = 'OVERFLOW - NO COMPUTATION, ';
var DEFAULT_OVERFLOW = 'FNU IS TOO LARGE OR CABS(Z) IS TOO SMALL OR BOTH';
var HALF_ACCURACY = 'CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE BUT LOSSES OF' + ' SIGNIFICANCE BY ARGUMENT REDUCTION PRODUCE LESS THAN HALF OF MACHINE ACCURACY.';
var TOO_LARGE = 'CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTATION BECAUSE OF' + ' COMPLETE LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION.';
var TERMINATION = 'ERROR - NO COMPUTATION ALGORITHM TERMINATION CONDITION NOT MET';

function besselCommon(nu, z, func, overflowMessage) {
 var hankelType = arguments.length > 4 && arguments[4] !== undefined ? arguments[4] : null;

 var nz = void 0,
 ierr = void 0;
 z = _complex2.default.ensureComplex(z);
 var n = 1;
 var cyr = new Array(1);
 var cyi = new Array(1);
 if (hankelType) {
 var _func = func(z[0], z[1], nu, 1, hankelType, n, cyr, cyi);

 var _func2 = _slicedToArray(_func, 2);

 nz = _func2[0];
 ierr = _func2[1];
 } else {
 var _func3 = func(z[0], z[1], nu, 1, n, cyr, cyi);

 var _func4 = _slicedToArray(_func3, 2);

 nz = _func4[0];
 ierr = _func4[1];
 }

 checkUnderflow(nz, func);
 checkErrorCode(ierr, overflowMessage, func);
 return [cyr[0], cyi[0]];
}

function airyCommon(z, func, overflowMessage) {
 z = _complex2.default.ensureComplex(z);

 var _func5 = func(z[0], z[1], 0, 1),
 _func6 = _slicedToArray(_func5, 4),
 air = _func6[0],
 aii = _func6[1],
 nz = _func6[2],
 ierr = _func6[3];

 checkUnderflow(nz, func);
 checkErrorCode(ierr, overflowMessage, func);
 return [air, aii];
}

function besselJ(nu, z) {
 var overflowMessage = 'AIMAG(Z) TOO LARGE ON KODE=1';
 return besselCommon(nu, z, _zbesj.zbesj, overflowMessage);
}

function besselY(nu, z) {
 var overflowMessage = DEFAULT_OVERFLOW;
 return besselCommon(nu, z, _zbesy.zbesy, overflowMessage);
}

function besselI(nu, z) {
 var overflowMessage = 'REAL(Z) TOO LARGE ON KODE=1';
 return besselCommon(nu, z, _zbesi.zbesi, overflowMessage);
}

function besselK(nu, z) {
 var overflowMessage = DEFAULT_OVERFLOW;
 return besselCommon(nu, z, _zbesk.zbesk, overflowMessage);
}

// Hankel functions type 1 or type 2
function besselH(nu, z, type) {
 var overflowMessage = DEFAULT_OVERFLOW;
 return besselCommon(nu, z, _zbesh.zbesh, overflowMessage, type);
}

function airyA(z) {
 var overflowMessage = 'REAL(ZTA) TOO LARGE ON KODE=1';
 return airyCommon(z, _zairy.zairy, overflowMessage);
}

function airyB(z) {
 var overflowMessage = 'REAL(Z) TOO LARGE ON KODE=1';
 return airyCommon(z, _zbiry.zbiry, overflowMessage);
}

function checkUnderflow(nz, func) {
 if (nz !== 0) {
 // console.warn(func.name + ' returned nonzero number of underflows');
 }
}

function checkErrorCode(ierr, overflowMessage, func) {
 switch (ierr) {
 case 0:
 // NORMAL RETURN - COMPUTATION COMPLETED
 break;
 case 1:
 throw new Error(INPUT_ERROR);
 case 2:
 throw new Error(OVERFLOW + overflowMessage);
 case 3:
 throw new Error(HALF_ACCURACY);
 case 4:
 throw new Error(TOO_LARGE);
 case 5:
 throw new Error(TERMINATION);
 default:
 throw new Error('Unexpected error code from ' + func.name);
 }
}

exports.besselJ = besselJ;
exports.besselY = besselY;
exports.besselI = besselI;
exports.besselK = besselK;
exports.besselH = besselH;
exports.airyA = airyA;
exports.airyB = airyB;
},{"../utils/complex.js":89,"./amos/zairy.js":14,"./amos/zbesh.js":16,"./amos/zbesi.js":17,"./amos/zbesj.js":18,"./amos/zbesk.js":19,"./amos/zbesy.js":20,"./amos/zbiry.js":22}],10:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.dgamln = dgamln;

var _d1mach = require('../../utils/fortran-utils/d1mach.js');

var _i1mach = require('../../utils/fortran-utils/i1mach.js');

// DOUBLE PRECISION FUNCTION DGAMLN(Z,IERR)
// ***BEGIN PROLOGUE DGAMLN
// ***DATE WRITTEN 830501 (YYMMDD)
// ***REVISION DATE 830501 (YYMMDD)
// ***PORT TO ECMASCRIPT 201801 (YYYYMM)
// ***CATEGORY NO. B5F
// ***KEYWORDS GAMMA FUNCTION,LOGARITHM OF GAMMA FUNCTION
// ***AUTHOR (FORTRAN) AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
// ***AUTHOR (ECMASCRIPT) ERB, KC, KINGS DISTRIBUTED SYSTEMS
// ***PURPOSE TO COMPUTE THE LOGARITHM OF THE GAMMA FUNCTION
// ***DESCRIPTION
//
// **** A DOUBLE PRECISION ROUTINE ****
// DGAMLN COMPUTES THE NATURAL LOG OF THE GAMMA FUNCTION FOR
// Z.GT.0. THE ASYMPTOTIC EXPANSION IS USED TO GENERATE VALUES
// GREATER THAN ZMIN WHICH ARE ADJUSTED BY THE RECURSION
// G(Z+1)=Z*G(Z) FOR Z.LE.ZMIN. THE FUNCTION WAS MADE AS
// PORTABLE AS POSSIBLE BY COMPUTIMG ZMIN FROM THE NUMBER OF BASE
// 10 DIGITS IN A WORD, RLN=AMAX1(-ALOG10(R1MACH(4)),0.5E-18)
// LIMITED TO 18 DIGITS OF (RELATIVE) ACCURACY.
//
// SINCE INTEGER ARGUMENTS ARE COMMON, A TABLE LOOK UP ON 100
// VALUES IS USED FOR SPEED OF EXECUTION.
//
// DESCRIPTION OF ARGUMENTS
//
// INPUT Z IS D0UBLE PRECISION
// Z - ARGUMENT, Z.GT.0.0D0
//
// OUTPUT DGAMLN IS DOUBLE PRECISION
// DGAMLN - NATURAL LOG OF THE GAMMA FUNCTION AT Z.NE.0.0D0
//
//
// ***REFERENCES COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
// BY D. E. AMOS, SAND83-0083, MAY, 1983.
// ***ROUTINES CALLED I1MACH,D1MACH
// ***END PROLOGUE DGAMLN
function dgamln(z) {
 var cf = void 0,
 con = void 0,
 fln = void 0,
 fz = void 0,
 gln = void 0,
 rln = void 0,
 s = void 0,
 tlg = void 0,
 trm = void 0,
 tst = void 0,
 t1 = void 0,
 wdtol = void 0,
 zdmy = void 0,
 zinc = void 0,
 zm = void 0,
 zmin = void 0,
 zp = void 0,
 zsq = void 0,
 i = void 0,
 i1m = void 0,
 k = void 0,
 mz = void 0,
 nz = void 0;
 gln = [0.00000000000000000e+00, 0.00000000000000000e+00, 6.93147180559945309e-01, 1.79175946922805500e+00, 3.17805383034794562e+00, 4.78749174278204599e+00, 6.57925121201010100e+00, 8.52516136106541430e+00, 1.06046029027452502e+01, 1.28018274800814696e+01, 1.51044125730755153e+01, 1.75023078458738858e+01, 1.99872144956618861e+01, 2.25521638531234229e+01, 2.51912211827386815e+01, 2.78992713838408916e+01, 3.06718601060806728e+01, 3.35050734501368889e+01, 3.63954452080330536e+01, 3.93398841871994940e+01, 4.23356164607534850e+01, 4.53801388984769080e+01, 4.84711813518352239e+01, 5.16066755677643736e+01, 5.47847293981123192e+01, 5.80036052229805199e+01, 6.12617017610020020e+01, 6.45575386270063311e+01, 6.78897431371815350e+01, 7.12570389671680090e+01, 7.46582363488301644e+01, 7.80922235533153106e+01, 8.15579594561150372e+01, 8.50544670175815174e+01, 8.85808275421976788e+01, 9.21361756036870925e+01, 9.57196945421432025e+01, 9.93306124547874269e+01, 1.02968198614513813e+02, 1.06631760260643459e+02, 1.10320639714757395e+02, 1.14034211781461703e+02, 1.17771881399745072e+02, 1.21533081515438634e+02, 1.25317271149356895e+02, 1.29123933639127215e+02, 1.32952575035616310e+02, 1.36802722637326368e+02, 1.40673923648234259e+02, 1.44565743946344886e+02, 1.48477766951773032e+02, 1.52409592584497358e+02, 1.56360836303078785e+02, 1.60331128216630907e+02, 1.64320112263195181e+02, 1.68327445448427652e+02, 1.72352797139162802e+02, 1.76395848406997352e+02, 1.80456291417543771e+02, 1.84533828861449491e+02, 1.88628173423671591e+02, 1.92739047287844902e+02, 1.96866181672889994e+02, 2.01009316399281527e+02, 2.05168199482641199e+02, 2.09342586752536836e+02, 2.13532241494563261e+02, 2.17736934113954227e+02, 2.21956441819130334e+02, 2.26190548323727593e+02, 2.30439043565776952e+02, 2.34701723442818268e+02, 2.38978389561834323e+02, 2.43268849002982714e+02, 2.47572914096186884e+02, 2.51890402209723194e+02, 2.56221135550009525e+02, 2.60564940971863209e+02, 2.64921649798552801e+02, 2.69291097651019823e+02, 2.73673124285693704e+02, 2.78067573440366143e+02, 2.82474292687630396e+02, 2.86893133295426994e+02, 2.91323950094270308e+02, 2.95766601350760624e+02, 3.00220948647014132e+02, 3.04686856765668715e+02, 3.09164193580146922e+02, 3.13652829949879062e+02, 3.18152639620209327e+02, 3.22663499126726177e+02, 3.27185287703775217e+02, 3.31717887196928473e+02, 3.36261181979198477e+02, 3.40815058870799018e+02, 3.45379407062266854e+02, 3.49954118040770237e+02, 3.54539085519440809e+02, 3.59134205369575399e+02];
 cf = [8.33333333333333333e-02, -2.77777777777777778e-03, 7.93650793650793651e-04, -5.95238095238095238e-04, 8.41750841750841751e-04, -1.91752691752691753e-03, 6.41025641025641026e-03, -2.95506535947712418e-02, 1.79644372368830573e-01, -1.39243221690590112e+00, 1.34028640441683920e+01, -1.56848284626002017e+02, 2.19310333333333333e+03, -3.61087712537249894e+04, 6.91472268851313067e+05, -1.52382215394074162e+07, 3.82900751391414141e+08, -1.08822660357843911e+10, 3.47320283765002252e+11, -1.23696021422692745e+13, 4.88788064793079335e+14, -2.13203339609193739e+16];
 con = 1.83787706640934548;
 // c***first executable statement dgamln
 // ierr = 0;
 if (z <= 0.0) {
 // ierr = 1;
 return null;
 }
 if (z > 101.0) {
 // go to 10
 } else {
 nz = Math.trunc(z);
 fz = z - nz;
 if (fz > 0.0 || nz > 100) {
 //
 } else {
 return gln[nz - 1];
 }
 }
 // 10 continue
 wdtol = (0, _d1mach.d1mach)(4);
 wdtol = Math.max(wdtol, 0.5e-18);
 i1m = (0, _i1mach.i1mach)(14);
 rln = (0, _d1mach.d1mach)(5) * i1m;
 fln = Math.min(rln, 20.0);
 fln = Math.max(fln, 3.0);
 fln = fln - 3.0;
 zm = 1.8000 + 0.3875 * fln;
 mz = Math.trunc(zm) + 1;
 zmin = mz;
 zdmy = z;
 zinc = 0.0;
 if (z >= zmin) {
 // go to 20
 } else {
 zinc = zmin - nz;
 zdmy = z + zinc;
 }
 // 20 continue
 zp = 1.0 / zdmy;
 t1 = cf[0] * zp;
 s = t1;
 if (zp < wdtol) {
 // go to 40
 } else {
 zsq = zp * zp;
 tst = t1 * wdtol;
 // do 30 k=2,22
 for (k = 1; k < 22; k++) {
 zp = zp * zsq;
 trm = cf[k] * zp;
 if (Math.abs(trm) < tst) break; // go to 40
 s = s + trm;
 }
 // 30 continue
 }
 // 40 continue
 if (zinc !== 0.0) {
 // go to 50
 } else {
 tlg = Math.log(z);
 return z * (tlg - 1.0) + 0.5 * (con - tlg) + s;
 }
 // 50 continue
 zp = 1.0;
 nz = Math.trunc(zinc);
 // do 60 i=1,nz
 for (i = 1; i <= nz; i++) {
 zp = zp * (z + (i - 1));
 }
 // 60 continue
 tlg = Math.log(zdmy);
 return zdmy * (tlg - 1.0) - Math.log(zp) + 0.5 * (con - tlg) + s;
}
},{"../../utils/fortran-utils/d1mach.js":91,"../../utils/fortran-utils/i1mach.js":92}],11:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.azabs = azabs;
// DOUBLE PRECISION FUNCTION AZABS(ZR, ZI)
// ***BEGIN PROLOGUE AZABS
// ***REFER TO ZBESH,ZBESI,ZBESJ,ZBESK,ZBESY,ZAIRY,ZBIRY
//
// AZABS COMPUTES THE ABSOLUTE VALUE OR MAGNITUDE OF A DOUBLE
// PRECISION COMPLEX VARIABLE CMPLX(ZR,ZI)
//
// ***ROUTINES CALLED (NONE)
// ***END PROLOGUE AZABS
function azabs(zr, zi) {
 var u = void 0,
 v = void 0,
 q = void 0,
 s = void 0;
 u = Math.abs(zr);
 v = Math.abs(zi);
 s = u + v;

 if (s === 0) return 0;
 if (u > v) {
 q = v / u;
 return u * Math.sqrt(1 + q * q);
 }
 q = u / v;
 return v * Math.sqrt(1 + q * q);
}
},{}],12:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});

var _slicedToArray = function () { function sliceIterator(arr, i) { var _arr = []; var _n = true; var _d = false; var _e = undefined; try { for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) { _arr.push(_s.value); if (i && _arr.length === i) break; } } catch (err) { _d = true; _e = err; } finally { try { if (!_n && _i["return"]) _i["return"](); } finally { if (_d) throw _e; } } return _arr; } return function (arr, i) { if (Array.isArray(arr)) { return arr; } else if (Symbol.iterator in Object(arr)) { return sliceIterator(arr, i); } else { throw new TypeError("Invalid attempt to destructure non-iterable instance"); } }; }(); /* eslint-disable no-labels */
/* eslint-disable no-fallthrough */
// SUBROUTINE ZACAI(ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, RL, TOL,
// * ELIM, ALIM)
// ***BEGIN PROLOGUE ZACAI
// ***REFER TO ZAIRY
//
// ZACAI APPLIES THE ANALYTIC CONTINUATION FORMULA
//
// K(FNU,ZN*EXP(MP))=K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN)
// MP=PI*MR*CMPLX(0.0,1.0)
//
// TO CONTINUE THE K FUNCTION FROM THE RIGHT HALF TO THE LEFT
// HALF Z PLANE FOR USE WITH ZAIRY WHERE FNU=1/3 OR 2/3 AND N=1.
// ZACAI IS THE SAME AS ZACON WITH THE PARTS FOR LARGER ORDERS AND
// RECURRENCE REMOVED. A RECURSIVE CALL TO ZACON CAN RESULT IF ZACON
// IS CALLED FROM ZAIRY.
//
// ***ROUTINES CALLED ZASYI,ZBKNU,ZMLRI,ZSERI,ZS1S2,D1MACH,AZABS
// ***END PROLOGUE ZACAI


exports.zacai = zacai;

var _zasyi = require('./zasyi.js');

var _zbknu = require('./zbknu.js');

var _zmlri = require('./zmlri.js');

var _zseri = require('./zseri.js');

var _zs1s3 = require('./zs1s2.js');

var _d1mach = require('../../utils/fortran-utils/d1mach.js');

var _zabs = require('./zabs.js');

var _fortranHelpers = require('../../utils/fortranHelpers.js');

var ft = _interopRequireWildcard(_fortranHelpers);

function _interopRequireWildcard(obj) { if (obj && obj.__esModule) { return obj; } else { var newObj = {}; if (obj != null) { for (var key in obj) { if (Object.prototype.hasOwnProperty.call(obj, key)) newObj[key] = obj[key]; } } newObj.default = obj; return newObj; } }

function zacai(zr, zi, fnu, kode, mr, n, yr, yi, rl, tol, elim, alim) {
 var arg = void 0,
 ascle = void 0,
 az = void 0,
 csgnr = void 0,
 csgni = void 0,
 cspnr = void 0,
 cspni = void 0,
 c1r = void 0,
 c1i = void 0,
 c2r = void 0,
 c2i = void 0,
 cyr = void 0,
 cyi = void 0,
 dfnu = void 0,
 fmr = void 0,
 pi = void 0,
 sgn = void 0,
 yy = void 0,
 znr = void 0,
 zni = void 0,
 inu = void 0,
 iuf = void 0,
 nn = void 0,
 nw = void 0,
 nz = void 0;
 cyr = new Array(2);
 cyi = new Array(2);
 pi = 3.14159265358979324;

 var goToLabel = 0;
 mainExecutionLoop: while (true) {
 switch (goToLabel) {
 case 0:
 nz = 0;
 znr = -zr;
 zni = -zi;
 az = (0, _zabs.azabs)(zr, zi);
 nn = n;
 dfnu = fnu + (n - 1);
 if (az <= 2.0) {
 goToLabel = 10;break;
 }
 if (az * az * 0.25 > dfnu + 1.0) {
 goToLabel = 20;break;
 }
 case 10:
 // c-----------------------------------------------------------------------
 // c power series for the i function
 // c-----------------------------------------------------------------------
 nw = (0, _zseri.zseri)(znr, zni, fnu, kode, nn, yr, yi, tol, elim, alim);
 goToLabel = 40;break;
 case 20:
 if (az < rl) {
 goToLabel = 30;break;
 }
 // c-----------------------------------------------------------------------
 // c asymptotic expansion for large z for the i function
 // c-----------------------------------------------------------------------
 nw = (0, _zasyi.zasyi)(znr, zni, fnu, kode, nn, yr, yi, rl, tol, elim, alim);
 if (nw < 0) {
 goToLabel = 80;break;
 }
 goToLabel = 40;break;
 case 30:
 // c-----------------------------------------------------------------------
 // c miller algorithm normalized by the series for the i function
 // c-----------------------------------------------------------------------
 nw = (0, _zmlri.zmlri)(znr, zni, fnu, kode, nn, yr, yi, tol);
 if (nw < 0) {
 goToLabel = 80;break;
 }
 case 40:
 // c-----------------------------------------------------------------------
 // c analytic continuation to the left half plane for the k function
 // c-----------------------------------------------------------------------
 nw = (0, _zbknu.zbknu)(znr, zni, fnu, kode, 1, cyr, cyi, tol, elim, alim);
 if (nw !== 0) {
 goToLabel = 80;break;
 }
 fmr = mr;
 sgn = -ft.sign(pi, fmr);
 csgnr = 0.0;
 csgni = sgn;
 if (kode === 1) {
 goToLabel = 50;break;
 }
 yy = -zni;
 csgnr = -csgni * Math.sin(yy);
 csgni = csgni * Math.cos(yy);
 case 50:
 // c-----------------------------------------------------------------------
 // c calculate cspn=exp(fnu*pi*i) to minimize losses of significance
 // c when fnu is large
 // c-----------------------------------------------------------------------
 inu = Math.trunc(fnu);
 arg = (fnu - inu) * sgn;
 cspnr = Math.cos(arg);
 cspni = Math.sin(arg);
 if (inu % 2 === 0) {
 goToLabel = 60;break;
 }
 cspnr = -cspnr;
 cspni = -cspni;
 case 60:
 c1r = cyr[0];
 c1i = cyi[0];
 c2r = yr[0];
 c2i = yi[0];
 if (kode === 1) {
 goToLabel = 70;break;
 }
 iuf = 0;
 ascle = 1.0e+3 * (0, _d1mach.d1mach)(1) / tol;

 var _zs1s = (0, _zs1s3.zs1s2)(znr, zni, c1r, c1i, c2r, c2i, ascle, alim, iuf);

 var _zs1s2 = _slicedToArray(_zs1s, 6);

 c1r = _zs1s2[0];
 c1i = _zs1s2[1];
 c2r = _zs1s2[2];
 c2i = _zs1s2[3];
 nw = _zs1s2[4];
 iuf = _zs1s2[5];

 nz = nz + nw;
 case 70:
 yr[0] = cspnr * c1r - cspni * c1i + csgnr * c2r - csgni * c2i;
 yi[0] = cspnr * c1i + cspni * c1r + csgnr * c2i + csgni * c2r;
 break mainExecutionLoop;
 case 80:
 nz = -1;
 if (nw === -2) nz = -2;
 default:
 break mainExecutionLoop;
 }
 }

 return nz;
}
},{"../../utils/fortran-utils/d1mach.js":91,"../../utils/fortranHelpers.js":93,"./zabs.js":11,"./zasyi.js":15,"./zbknu.js":23,"./zmlri.js":30,"./zs1s2.js":33,"./zseri.js":34}],13:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});

var _slicedToArray = function () { function sliceIterator(arr, i) { var _arr = []; var _n = true; var _d = false; var _e = undefined; try { for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) { _arr.push(_s.value); if (i && _arr.length === i) break; } } catch (err) { _d = true; _e = err; } finally { try { if (!_n && _i["return"]) _i["return"](); } finally { if (_d) throw _e; } } return _arr; } return function (arr, i) { if (Array.isArray(arr)) { return arr; } else if (Symbol.iterator in Object(arr)) { return sliceIterator(arr, i); } else { throw new TypeError("Invalid attempt to destructure non-iterable instance"); } }; }(); /* eslint-disable no-labels */
/* eslint-disable no-fallthrough */
// SUBROUTINE ZACON(ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, RL, FNUL,
// * TOL, ELIM, ALIM)
// ***BEGIN PROLOGUE ZACON
// ***REFER TO ZBESK,ZBESH
//
// ZACON APPLIES THE ANALYTIC CONTINUATION FORMULA
//
// K(FNU,ZN*EXP(MP))=K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN)
// MP=PI*MR*CMPLX(0.0,1.0)
//
// TO CONTINUE THE K FUNCTION FROM THE RIGHT HALF TO THE LEFT
// HALF Z PLANE
//
// ***ROUTINES CALLED ZBINU,ZBKNU,ZS1S2,D1MACH,AZABS,ZMLT
// ***END PROLOGUE ZACON
// COMPLEX CK,CONE,CSCL,CSCR,CSGN,CSPN,CY,CZERO,C1,C2,RZ,SC1,SC2,ST,
// *S1,S2,Y,Z,ZN


exports.zacon = zacon;

var _zbknu = require('./zbknu.js');

var _zbinu = require('./zbinu.js');

var _zs1s7 = require('./zs1s2.js');

var _d1mach = require('../../utils/fortran-utils/d1mach.js');

var _zabs = require('./zabs.js');

var _zmlt11 = require('./zmlt.js');

var _fortranHelpers = require('../../utils/fortranHelpers.js');

var ft = _interopRequireWildcard(_fortranHelpers);

function _interopRequireWildcard(obj) { if (obj && obj.__esModule) { return obj; } else { var newObj = {}; if (obj != null) { for (var key in obj) { if (Object.prototype.hasOwnProperty.call(obj, key)) newObj[key] = obj[key]; } } newObj.default = obj; return newObj; } }

function zacon(zr, zi, fnu, kode, mr, n, yr, yi, rl, fnul, tol, elim, alim) {
 var arg = void 0,
 ascle = void 0,
 as2 = void 0,
 azn = void 0,
 bry = void 0,
 bscle = void 0,
 cki = void 0,
 ckr = void 0,
 coner = void 0,
 cpn = void 0,
 cscl = void 0,
 cscr = void 0,
 csgni = void 0,
 csgnr = void 0,
 cspni = void 0,
 cspnr = void 0,
 csr = void 0,
 csrr = void 0,
 cssr = void 0,
 cyi = void 0,
 cyr = void 0,
 c1i = void 0,
 c1m = void 0,
 c1r = void 0,
 c2i = void 0,
 c2r = void 0,
 fmr = void 0,
 fn = void 0,
 pi = void 0,
 pti = void 0,
 ptr = void 0,
 razn = void 0,
 rzi = void 0,
 rzr = void 0,
 sc1i = void 0,
 sc1r = void 0,
 sc2i = void 0,
 sc2r = void 0,
 sgn = void 0,
 spn = void 0,
 sti = void 0,
 str = void 0,
 s1i = void 0,
 s1r = void 0,
 s2i = void 0,
 s2r = void 0,
 yy = void 0,
 zeror = void 0,
 zni = void 0,
 znr = void 0,
 i = void 0,
 inu = void 0,
 iuf = void 0,
 kflag = void 0,
 nn = void 0,
 nw = void 0,
 nz = void 0;

 cyr = new Float64Array(2);
 cyi = new Float64Array(2);
 cssr = new Float64Array(3);
 csrr = new Float64Array(3);
 bry = new Float64Array(3);
 pi = 3.14159265358979324;

 zeror = 0;
 coner = 1;
 var goToLabel = 0;
 mainExecutionLoop: while (true) {
 switch (goToLabel) {
 case 0:
 nz = 0;
 znr = -zr;
 zni = -zi;
 nn = n;
 nw = (0, _zbinu.zbinu)(znr, zni, fnu, kode, nn, yr, yi, rl, fnul, tol, elim, alim);
 if (nw < 0) {
 goToLabel = 90;break;
 }
 // c-----------------------------------------------------------------------
 // c analytic continuation to the left half plane for the k function
 // c-----------------------------------------------------------------------
 nn = Math.min(2, n);
 nw = (0, _zbknu.zbknu)(znr, zni, fnu, kode, nn, cyr, cyi, tol, elim, alim);
 if (nw !== 0) {
 goToLabel = 90;break;
 }
 s1r = cyr[0];
 s1i = cyi[0];
 fmr = mr;
 sgn = -ft.sign(pi, fmr);
 csgnr = zeror;
 csgni = sgn;
 if (kode === 1) {
 goToLabel = 10;break;
 }
 yy = -zni;
 cpn = Math.cos(yy);
 spn = Math.sin(yy);

 var _zmlt = (0, _zmlt11.zmlt)(csgnr, csgni, cpn, spn);

 var _zmlt2 = _slicedToArray(_zmlt, 2);

 csgnr = _zmlt2[0];
 csgni = _zmlt2[1];

 case 10:
 // c-----------------------------------------------------------------------
 // c calculate cspn=exp(fnu*pi*i) to minimize losses of significance
 // c when fnu is large
 // c-----------------------------------------------------------------------
 inu = Math.trunc(fnu);
 arg = (fnu - inu) * sgn;
 cpn = Math.cos(arg);
 spn = Math.sin(arg);
 cspnr = cpn;
 cspni = spn;
 if (inu % 2 === 0) {
 goToLabel = 20;break;
 }
 cspnr = -cspnr;
 cspni = -cspni;
 case 20:
 iuf = 0;
 c1r = s1r;
 c1i = s1i;
 c2r = yr[0];
 c2i = yi[0];
 ascle = 1.0e+3 * (0, _d1mach.d1mach)(1) / tol;
 if (kode === 1) {
 goToLabel = 30;break;
 }

 var _zs1s = (0, _zs1s7.zs1s2)(znr, zni, c1r, c1i, c2r, c2i, ascle, alim, iuf);

 var _zs1s2 = _slicedToArray(_zs1s, 6);

 nw = _zs1s2[0];
 c1r = _zs1s2[1];
 c1i = _zs1s2[2];
 c2r = _zs1s2[3];
 c2i = _zs1s2[4];
 iuf = _zs1s2[5];

 nz = nz + nw;
 sc1r = c1r;
 sc1i = c1i;
 case 30:
 var _zmlt3 = (0, _zmlt11.zmlt)(cspnr, cspni, c1r, c1i);

 var _zmlt4 = _slicedToArray(_zmlt3, 2);

 str = _zmlt4[0];
 sti = _zmlt4[1];

 var _zmlt5 = (0, _zmlt11.zmlt)(csgnr, csgni, c2r, c2i);

 var _zmlt6 = _slicedToArray(_zmlt5, 2);

 ptr = _zmlt6[0];
 pti = _zmlt6[1];

 yr[0] = str + ptr;
 yi[0] = sti + pti;
 if (n === 1) break mainExecutionLoop;
 cspnr = -cspnr;
 cspni = -cspni;
 s2r = cyr[1];
 s2i = cyi[1];
 c1r = s2r;
 c1i = s2i;
 c2r = yr[1];
 c2i = yi[1];
 if (kode === 1) {
 goToLabel = 40;break;
 }

 var _zs1s3 = (0, _zs1s7.zs1s2)(znr, zni, c1r, c1i, c2r, c2i, ascle, alim, iuf);

 var _zs1s4 = _slicedToArray(_zs1s3, 6);

 nw = _zs1s4[0];
 c1r = _zs1s4[1];
 c1i = _zs1s4[2];
 c2r = _zs1s4[3];
 c2i = _zs1s4[4];
 iuf = _zs1s4[5];

 nz = nz + nw;
 sc2r = c1r;
 sc2i = c1i;
 case 40:
 var _zmlt7 = (0, _zmlt11.zmlt)(cspnr, cspni, c1r, c1i);

 var _zmlt8 = _slicedToArray(_zmlt7, 2);

 str = _zmlt8[0];
 sti = _zmlt8[1];

 var _zmlt9 = (0, _zmlt11.zmlt)(csgnr, csgni, c2r, c2i);

 var _zmlt10 = _slicedToArray(_zmlt9, 2);

 ptr = _zmlt10[0];
 pti = _zmlt10[1];

 yr[1] = str + ptr;
 yi[1] = sti + pti;
 if (n === 2) break mainExecutionLoop;
 cspnr = -cspnr;
 cspni = -cspni;
 azn = (0, _zabs.azabs)(znr, zni);
 razn = 1.0 / azn;
 str = znr * razn;
 sti = -zni * razn;
 rzr = (str + str) * razn;
 rzi = (sti + sti) * razn;
 fn = fnu + 1.0;
 ckr = fn * rzr;
 cki = fn * rzi;
 // c-----------------------------------------------------------------------
 // c scale near exponent extremes during recurrence on k functions
 // c-----------------------------------------------------------------------
 cscl = 1.0 / tol;
 cscr = tol;
 cssr[0] = cscl;
 cssr[1] = coner;
 cssr[2] = cscr;
 csrr[0] = cscr;
 csrr[1] = coner;
 csrr[2] = cscl;
 bry[0] = ascle;
 bry[1] = 1.0 / ascle;
 bry[2] = (0, _d1mach.d1mach)(2);
 as2 = (0, _zabs.azabs)(s2r, s2i);
 kflag = 2;
 if (as2 > bry[0]) {
 goToLabel = 50;break;
 }
 kflag = 1;
 goToLabel = 60;break;
 case 50:
 if (as2 < bry[1]) {
 goToLabel = 60;break;
 }
 kflag = 3;
 case 60:
 bscle = bry[kflag - 1];
 s1r = s1r * cssr[kflag - 1];
 s1i = s1i * cssr[kflag - 1];
 s2r = s2r * cssr[kflag - 1];
 s2i = s2i * cssr[kflag - 1];
 csr = csrr[kflag - 1];
 // do 80 i=3,n
 for (i = 3; i <= n; i++) {
 str = s2r;
 sti = s2i;
 s2r = ckr * str - cki * sti + s1r;
 s2i = ckr * sti + cki * str + s1i;
 s1r = str;
 s1i = sti;
 c1r = s2r * csr;
 c1i = s2i * csr;
 str = c1r;
 sti = c1i;
 c2r = yr[i - 1];
 c2i = yi[i - 1];
 if (kode === 1) {
 // go to 70
 } else {
 if (iuf < 0) {
 // go to 70
 } else {
 var _zs1s5 = (0, _zs1s7.zs1s2)(znr, zni, c1r, c1i, c2r, c2i, ascle, alim, iuf);

 var _zs1s6 = _slicedToArray(_zs1s5, 6);

 nw = _zs1s6[0];
 c1r = _zs1s6[1];
 c1i = _zs1s6[2];
 c2r = _zs1s6[3];
 c2i = _zs1s6[4];
 iuf = _zs1s6[5];

 nz = nz + nw;
 sc1r = sc2r;
 sc1i = sc2i;
 sc2r = c1r;
 sc2i = c1i;
 if (iuf !== 3) {
 // go to 70
 } else {
 iuf = -4;
 s1r = sc1r * cssr[kflag - 1];
 s1i = sc1i * cssr[kflag - 1];
 s2r = sc2r * cssr[kflag - 1];
 s2i = sc2i * cssr[kflag - 1];
 str = sc2r;
 sti = sc2i;
 }
 }
 }
 // 70 continue
 ptr = cspnr * c1r - cspni * c1i;
 pti = cspnr * c1i + cspni * c1r;
 yr[i - 1] = ptr + csgnr * c2r - csgni * c2i;
 yi[i - 1] = pti + csgnr * c2i + csgni * c2r;
 ckr = ckr + rzr;
 cki = cki + rzi;
 cspnr = -cspnr;
 cspni = -cspni;
 if (kflag >= 3) continue;
 ptr = Math.abs(c1r);
 pti = Math.abs(c1i);
 c1m = Math.max(ptr, pti);
 if (c1m <= bscle) continue;
 kflag = kflag + 1;
 bscle = bry[kflag - 1];
 s1r = s1r * csr;
 s1i = s1i * csr;
 s2r = str;
 s2i = sti;
 s1r = s1r * cssr[kflag - 1];
 s1i = s1i * cssr[kflag - 1];
 s2r = s2r * cssr[kflag - 1];
 s2i = s2i * cssr[kflag - 1];
 csr = csrr[kflag - 1];
 } // 80 continue
 break mainExecutionLoop;
 case 90:
 nz = -1;
 if (nw === -2) nz = -2;
 default:
 break mainExecutionLoop;
 }
 }
 return nz;
}
},{"../../utils/fortran-utils/d1mach.js":91,"../../utils/fortranHelpers.js":93,"./zabs.js":11,"./zbinu.js":21,"./zbknu.js":23,"./zmlt.js":31,"./zs1s2.js":33}],14:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});

var _slicedToArray = function () { function sliceIterator(arr, i) { var _arr = []; var _n = true; var _d = false; var _e = undefined; try { for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) { _arr.push(_s.value); if (i && _arr.length === i) break; } } catch (err) { _d = true; _e = err; } finally { try { if (!_n && _i["return"]) _i["return"](); } finally { if (_d) throw _e; } } return _arr; } return function (arr, i) { if (Array.isArray(arr)) { return arr; } else if (Symbol.iterator in Object(arr)) { return sliceIterator(arr, i); } else { throw new TypeError("Invalid attempt to destructure non-iterable instance"); } }; }(); /* eslint-disable no-labels */
/* eslint-disable no-fallthrough */
// SUBROUTINE ZAIRY(ZR, ZI, ID, KODE, AIR, AII, NZ, IERR)
// ***BEGIN PROLOGUE ZAIRY
// ***DATE WRITTEN 830501 (YYMMDD)
// ***REVISION DATE 890801 (YYMMDD)
// ***PORT TO ECMASCRIPT 201801 (YYYYMM)
// ***CATEGORY NO. B5K
// ***KEYWORDS AIRY FUNCTION,BESSEL FUNCTIONS OF ORDER ONE THIRD
// ***AUTHOR (FORTRAN) AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
// ***AUTHOR (ECMASCRIPT) ERB, KC, KINGS DISTRIBUTED SYSTEMS
// ***PURPOSE TO COMPUTE AIRY FUNCTIONS AI(Z) AND DAI(Z) FOR COMPLEX Z
// ***DESCRIPTION
//
// ***A DOUBLE PRECISION ROUTINE***
// ON KODE=1, ZAIRY COMPUTES THE COMPLEX AIRY FUNCTION AI(Z) OR
// ITS DERIVATIVE DAI(Z)/DZ ON ID=0 OR ID=1 RESPECTIVELY. ON
// KODE=2, A SCALING OPTION CEXP(ZTA)*AI(Z) OR CEXP(ZTA)*
// DAI(Z)/DZ IS PROVIDED TO REMOVE THE EXPONENTIAL DECAY IN
// -PI/3.LT.ARG(Z).LT.PI/3 AND THE EXPONENTIAL GROWTH IN
// PI/3.LT.ABS(ARG(Z)).LT.PI WHERE ZTA=(2/3)*Z*CSQRT(Z).
//
// WHILE THE AIRY FUNCTIONS AI(Z) AND DAI(Z)/DZ ARE ANALYTIC IN
// THE WHOLE Z PLANE, THE CORRESPONDING SCALED FUNCTIONS DEFINED
// FOR KODE=2 HAVE A CUT ALONG THE NEGATIVE REAL AXIS.
// DEFINTIONS AND NOTATION ARE FOUND IN THE NBS HANDBOOK OF
// MATHEMATICAL FUNCTIONS (REF. 1).
//
// INPUT ZR,ZI ARE DOUBLE PRECISION
// ZR,ZI - Z=CMPLX(ZR,ZI)
// ID - ORDER OF DERIVATIVE, ID=0 OR ID=1
// KODE - A PARAMETER TO INDICATE THE SCALING OPTION
// KODE= 1 RETURNS
// AI=AI(Z) ON ID=0 OR
// AI=DAI(Z)/DZ ON ID=1
// = 2 RETURNS
// AI=CEXP(ZTA)*AI(Z) ON ID=0 OR
// AI=CEXP(ZTA)*DAI(Z)/DZ ON ID=1 WHERE
// ZTA=(2/3)*Z*CSQRT(Z)
//
// OUTPUT AIR,AII ARE DOUBLE PRECISION
// AIR,AII- COMPLEX ANSWER DEPENDING ON THE CHOICES FOR ID AND
// KODE
// NZ - UNDERFLOW INDICATOR
// NZ= 0 , NORMAL RETURN
// NZ= 1 , AI=CMPLX(0.0D0,0.0D0) DUE TO UNDERFLOW IN
// -PI/3.LT.ARG(Z).LT.PI/3 ON KODE=1
// IERR - ERROR FLAG
// IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
// IERR=1, INPUT ERROR - NO COMPUTATION
// IERR=2, OVERFLOW - NO COMPUTATION, REAL(ZTA)
// TOO LARGE ON KODE=1
// IERR=3, CABS(Z) LARGE - COMPUTATION COMPLETED
// LOSSES OF SIGNIFCANCE BY ARGUMENT REDUCTION
// PRODUCE LESS THAN HALF OF MACHINE ACCURACY
// IERR=4, CABS(Z) TOO LARGE - NO COMPUTATION
// COMPLETE LOSS OF ACCURACY BY ARGUMENT
// REDUCTION
// IERR=5, ERROR - NO COMPUTATION,
// ALGORITHM TERMINATION CONDITION NOT MET
//
// ***LONG DESCRIPTION
//
// AI AND DAI ARE COMPUTED FOR CABS(Z).GT.1.0 FROM THE K BESSEL
// FUNCTIONS BY
//
// AI(Z)=C*SQRT(Z)*K(1/3,ZTA) , DAI(Z)=-C*Z*K(2/3,ZTA)
// C=1.0/(PI*SQRT(3.0))
// ZTA=(2/3)*Z**(3/2)
//
// WITH THE POWER SERIES FOR CABS(Z).LE.1.0.
//
// IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
// MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z IS LARGE, LOSSES
// OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR. CONSEQUENTLY, IF
// THE MAGNITUDE OF ZETA=(2/3)*Z**1.5 EXCEEDS U1=SQRT(0.5/UR),
// THEN LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR
// FLAG IERR=3 IS TRIGGERED WHERE UR=DMAX1(D1MACH(4),1.0D-18) IS
// DOUBLE PRECISION UNIT ROUNDOFF LIMITED TO 18 DIGITS PRECISION.
// ALSO, IF THE MAGNITUDE OF ZETA IS LARGER THAN U2=0.5/UR, THEN
// ALL SIGNIFICANCE IS LOST AND IERR=4. IN ORDER TO USE THE INT
// FUNCTION, ZETA MUST BE FURTHER RESTRICTED NOT TO EXCEED THE
// LARGEST INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF ZETA
// MUST BE RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2,
// AND U3 ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE
// PRECISION ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE
// PRECISION ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMIT-
// ING IN THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT THE MAG-
// NITUDE OF Z CANNOT EXCEED 3.1E+4 IN SINGLE AND 2.1E+6 IN
// DOUBLE PRECISION ARITHMETIC. THIS ALSO MEANS THAT ONE CAN
// EXPECT TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES,
// NO DIGITS IN SINGLE PRECISION AND ONLY 7 DIGITS IN DOUBLE
// PRECISION ARITHMETIC. SIMILAR CONSIDERATIONS HOLD FOR OTHER
// MACHINES.
//
// THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
// BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
// ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
// SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
// ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
// ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
// CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
// HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
// ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
// SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
// THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
// 0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
// THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
// COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
// BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
// COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
// MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
// THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
// OR -PI/2+P.
//
// ***REFERENCES HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
// AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
// COMMERCE, 1955.
//
// COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
// AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
//
// A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
// ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
// 1018, MAY, 1985
//
// A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
// ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, TRANS.
// MATH. SOFTWARE, 1986
//
// ***ROUTINES CALLED ZACAI,ZBKNU,AZEXP,AZSQRT,I1MACH,D1MACH
// ***END PROLOGUE ZAIRY


exports.zairy = zairy;

var _zacai = require('./zacai.js');

var _zbknu = require('./zbknu.js');

var _zexp = require('./zexp.js');

var _zsqrt = require('./zsqrt.js');

var _zabs = require('./zabs.js');

var _i1mach = require('../../utils/fortran-utils/i1mach.js');

var _d1mach = require('../../utils/fortran-utils/d1mach.js');

function zairy(zr, zi, id, kode) {
 var aa = void 0,
 ad = void 0,
 aii = void 0,
 air = void 0,
 ak = void 0,
 alim = void 0,
 atrm = void 0,
 az = void 0,
 az3 = void 0,
 bk = void 0,
 cc = void 0,
 ck = void 0,
 coef = void 0,
 conei = void 0,
 coner = void 0,
 csqi = void 0,
 csqr = void 0,
 cyi = void 0,
 cyr = void 0,
 c1 = void 0,
 c2 = void 0,
 dig = void 0,
 dk = void 0,
 d1 = void 0,
 d2 = void 0,
 elim = void 0,
 fid = void 0,
 fnu = void 0,
 ptr = void 0,
 rl = void 0,
 r1m5 = void 0,
 sfac = void 0,
 sti = void 0,
 str = void 0,
 s1i = void 0,
 s1r = void 0,
 s2i = void 0,
 s2r = void 0,
 tol = void 0,
 trm1i = void 0,
 trm1r = void 0,
 trm2i = void 0,
 trm2r = void 0,
 tth = void 0,
 zeroi = void 0,
 zeror = void 0,
 ztai = void 0,
 ztar = void 0,
 z3i = void 0,
 z3r = void 0,
 alaz = void 0,
 bb = void 0,
 ierr = void 0,
 iflag = void 0,
 k = void 0,
 k1 = void 0,
 k2 = void 0,
 mr = void 0,
 nn = void 0,
 nz = void 0;
 cyr = [];
 cyi = [];
 tth = 6.66666666666666667e-01;
 c1 = 3.55028053887817240e-01;
 c2 = 2.58819403792806799e-01;
 coef = 1.83776298473930683e-01;
 zeror = 0;
 zeroi = 0;
 coner = 1;
 conei = 0;


 var goToLabel = 0;
 mainExecutionLoop: while (true) {
 switch (goToLabel) {
 case 0:
 ierr = 0;
 nz = 0;
 if (id < 0 || id > 1) ierr = 1;
 if (kode < 1 || kode > 2) ierr = 1;
 if (ierr !== 0) break mainExecutionLoop;
 az = (0, _zabs.azabs)(zr, zi);
 tol = Math.max((0, _d1mach.d1mach)(4), 1.0e-18);
 fid = id;
 if (az > 1.0) {
 goToLabel = 70;break;
 }
 // c-----------------------------------------------------------------------
 // c power series for cabs(z) <= 1.
 // c-----------------------------------------------------------------------
 s1r = coner;
 s1i = conei;
 s2r = coner;
 s2i = conei;
 if (az < tol) {
 goToLabel = 170;break;
 }
 aa = az * az;
 if (aa < tol / az) {
 goToLabel = 40;break;
 }
 trm1r = coner;
 trm1i = conei;
 trm2r = coner;
 trm2i = conei;
 atrm = 1.0;
 str = zr * zr - zi * zi;
 sti = zr * zi + zi * zr;
 z3r = str * zr - sti * zi;
 z3i = str * zi + sti * zr;
 az3 = az * aa;
 ak = 2.0 + fid;
 bk = 3.0 - fid - fid;
 ck = 4.0 - fid;
 dk = 3.0 + fid + fid;
 d1 = ak * dk;
 d2 = bk * ck;
 ad = Math.min(d1, d2);
 ak = 24.0 + 9.0 * fid;
 bk = 30.0 - 9.0 * fid;
 // do 30 k=1,25
 for (k = 1; k <= 25; k++) {
 str = (trm1r * z3r - trm1i * z3i) / d1;
 trm1i = (trm1r * z3i + trm1i * z3r) / d1;
 trm1r = str;
 s1r = s1r + trm1r;
 s1i = s1i + trm1i;
 str = (trm2r * z3r - trm2i * z3i) / d2;
 trm2i = (trm2r * z3i + trm2i * z3r) / d2;
 trm2r = str;
 s2r = s2r + trm2r;
 s2i = s2i + trm2i;
 atrm = atrm * az3 / ad;
 d1 = d1 + ak;
 d2 = d2 + bk;
 ad = Math.min(d1, d2);
 if (atrm < tol * ad) break; // go to 40
 ak = ak + 18.0;
 bk = bk + 18.0;
 }
 // 30 continue
 // 40 continue
 if (id === 1) {
 goToLabel = 50;break;
 }
 air = s1r * c1 - c2 * (zr * s2r - zi * s2i);
 aii = s1i * c1 - c2 * (zr * s2i + zi * s2r);
 if (kode === 1) break mainExecutionLoop;

 var _azsqrt = (0, _zsqrt.azsqrt)(zr, zi);

 var _azsqrt2 = _slicedToArray(_azsqrt, 2);

 str = _azsqrt2[0];
 sti = _azsqrt2[1];

 ztar = tth * (zr * str - zi * sti);
 ztai = tth * (zr * sti + zi * str);

 var _azexp = (0, _zexp.azexp)(ztar, ztai);

 var _azexp2 = _slicedToArray(_azexp, 2);

 str = _azexp2[0];
 sti = _azexp2[1];

 ptr = air * str - aii * sti;
 aii = air * sti + aii * str;
 air = ptr;
 break mainExecutionLoop;
 case 50:
 air = -s2r * c2;
 aii = -s2i * c2;
 if (az <= tol) {
 goToLabel = 60;break;
 }
 str = zr * s1r - zi * s1i;
 sti = zr * s1i + zi * s1r;
 cc = c1 / (1.0 + fid);
 air = air + cc * (str * zr - sti * zi);
 aii = aii + cc * (str * zi + sti * zr);
 case 60:
 if (kode === 1) break mainExecutionLoop;

 var _azsqrt3 = (0, _zsqrt.azsqrt)(zr, zi);

 var _azsqrt4 = _slicedToArray(_azsqrt3, 2);

 str = _azsqrt4[0];
 sti = _azsqrt4[1];

 ztar = tth * (zr * str - zi * sti);
 ztai = tth * (zr * sti + zi * str);

 var _azexp3 = (0, _zexp.azexp)(ztar, ztai);

 var _azexp4 = _slicedToArray(_azexp3, 2);

 str = _azexp4[0];
 sti = _azexp4[1];

 ptr = str * air - sti * aii;
 aii = str * aii + sti * air;
 air = ptr;
 break mainExecutionLoop;
 // c-----------------------------------------------------------------------
 // c case for cabs(z) > 1.0
 // c-----------------------------------------------------------------------
 case 70:
 fnu = (1.0 + fid) / 3.0;
 // c-----------------------------------------------------------------------
 // c set parameters related to machine constants.
 // c tol is the approximate unit roundoff limited to 1.0e-18.
 // c elim is the approximate exponential over- and underflow limit.
 // c exp(-elim) < exp(-alim)=exp(-elim)/tol and
 // c exp(elim) > exp(alim)=exp(elim)*tol are intervals near
 // c underflow and overflow limits where scaled arithmetic is done.
 // c rl is the lower boundary of the asymptotic expansion for large z.
 // c dig = number of base 10 digits in tol = 10**(-dig).
 // c-----------------------------------------------------------------------
 k1 = (0, _i1mach.i1mach)(15);
 k2 = (0, _i1mach.i1mach)(16);
 r1m5 = (0, _d1mach.d1mach)(5);
 k = Math.min(Math.abs(k1), Math.abs(k2));
 elim = 2.303 * (k * r1m5 - 3.0);
 k1 = (0, _i1mach.i1mach)(14) - 1;
 aa = r1m5 * k1;
 dig = Math.min(aa, 18.0);
 aa = aa * 2.303;
 alim = elim + Math.max(-aa, -41.45);
 rl = 1.2 * dig + 3.0;
 alaz = Math.log(az);
 // c--------------------------------------------------------------------------
 // c test for proper range
 // c-----------------------------------------------------------------------
 aa = 0.5 / tol;
 bb = (0, _i1mach.i1mach)(9) * 0.5;
 aa = Math.min(aa, bb);
 aa = aa ** tth;
 if (az > aa) {
 goToLabel = 260;break;
 }
 aa = Math.sqrt(aa);
 if (az > aa) ierr = 3;

 var _azsqrt5 = (0, _zsqrt.azsqrt)(zr, zi);

 var _azsqrt6 = _slicedToArray(_azsqrt5, 2);

 csqr = _azsqrt6[0];
 csqi = _azsqrt6[1];

 ztar = tth * (zr * csqr - zi * csqi);
 ztai = tth * (zr * csqi + zi * csqr);
 // c-----------------------------------------------------------------------
 // c re(zta) <= 0 when re(z) < 0, especially when im(z) is small
 // c-----------------------------------------------------------------------
 iflag = 0;
 sfac = 1.0;
 ak = ztai;
 if (zr >= 0.0) {
 goToLabel = 80;break;
 }
 bk = ztar;
 ck = -Math.abs(bk);
 ztar = ck;
 ztai = ak;
 case 80:
 if (zi !== 0.0) {
 goToLabel = 90;break;
 }
 if (zr > 0.0) {
 goToLabel = 90;break;
 }
 ztar = 0.0;
 ztai = ak;
 case 90:
 aa = ztar;
 if (aa >= 0.0 && zr > 0.0) {
 goToLabel = 110;break;
 }
 if (kode === 2) {
 goToLabel = 100;break;
 }
 // c-----------------------------------------------------------------------
 // c overflow test
 // c-----------------------------------------------------------------------
 if (aa > -alim) {
 goToLabel = 100;break;
 }
 aa = -aa + 0.25 * alaz;
 iflag = 1;
 sfac = tol;
 if (aa > elim) {
 goToLabel = 270;break;
 }
 case 100:
 // c-----------------------------------------------------------------------
 // c cbknu and cacon return exp(zta)*k(fnu,zta) on kode=2
 // c-----------------------------------------------------------------------
 mr = 1;
 if (zi < 0.0) mr = -1;
 nn = (0, _zacai.zacai)(ztar, ztai, fnu, kode, mr, 1, cyr, cyi, rl, tol, elim, alim);
 if (nn < 0) {
 goToLabel = 280;break;
 }
 nz = nz + nn;
 goToLabel = 130;break;
 case 110:
 if (kode === 2) {
 goToLabel = 120;break;
 }
 // c-----------------------------------------------------------------------
 // c underflow test
 // c-----------------------------------------------------------------------
 if (aa < alim) {
 goToLabel = 120;break;
 }
 aa = -aa - 0.25 * alaz;
 iflag = 2;
 sfac = 1.0 / tol;
 if (aa < -elim) {
 goToLabel = 210;break;
 }
 case 120:
 nz = (0, _zbknu.zbknu)(ztar, ztai, fnu, kode, 1, cyr, cyi, tol, elim, alim);
 case 130:
 s1r = cyr[0] * coef;
 s1i = cyi[0] * coef;
 if (iflag !== 0) {
 goToLabel = 150;break;
 }
 if (id === 1) {
 goToLabel = 140;break;
 }
 air = csqr * s1r - csqi * s1i;
 aii = csqr * s1i + csqi * s1r;
 break mainExecutionLoop;
 case 140:
 air = -(zr * s1r - zi * s1i);
 aii = -(zr * s1i + zi * s1r);
 break mainExecutionLoop;
 case 150:
 s1r = s1r * sfac;
 s1i = s1i * sfac;
 if (id === 1) {
 goToLabel = 160;break;
 }
 str = s1r * csqr - s1i * csqi;
 s1i = s1r * csqi + s1i * csqr;
 s1r = str;
 air = s1r / sfac;
 aii = s1i / sfac;
 break mainExecutionLoop;
 case 160:
 str = -(s1r * zr - s1i * zi);
 s1i = -(s1r * zi + s1i * zr);
 s1r = str;
 air = s1r / sfac;
 aii = s1i / sfac;
 break mainExecutionLoop;
 case 170:
 aa = 1.0e+3 * (0, _d1mach.d1mach)(1);
 s1r = zeror;
 s1i = zeroi;
 if (id === 1) {
 goToLabel = 190;break;
 }
 if (az <= aa) {
 goToLabel = 180;break;
 }
 s1r = c2 * zr;
 s1i = c2 * zi;
 case 180:
 air = c1 - s1r;
 aii = -s1i;
 break mainExecutionLoop;
 case 190:
 air = -c2;
 aii = 0.0;
 aa = Math.sqrt(aa);
 if (az <= aa) {
 goToLabel = 200;break;
 }
 s1r = 0.5 * (zr * zr - zi * zi);
 s1i = zr * zi;
 case 200:
 air = air + c1 * s1r;
 aii = aii + c1 * s1i;
 break mainExecutionLoop;
 case 210:
 nz = 1;
 air = zeror;
 aii = zeroi;
 break mainExecutionLoop;
 case 270:
 nz = 0;
 ierr = 2;
 break mainExecutionLoop;
 case 280:
 if (nn === -1) {
 goToLabel = 270;break;
 }
 nz = 0;
 ierr = 5;
 break mainExecutionLoop;
 case 260:
 ierr = 4;
 nz = 0;
 default:
 break mainExecutionLoop;
 }
 }

 return [air, aii, nz, ierr];
}
},{"../../utils/fortran-utils/d1mach.js":91,"../../utils/fortran-utils/i1mach.js":92,"./zabs.js":11,"./zacai.js":12,"./zbknu.js":23,"./zexp.js":27,"./zsqrt.js":36}],15:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});

var _slicedToArray = function () { function sliceIterator(arr, i) { var _arr = []; var _n = true; var _d = false; var _e = undefined; try { for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) { _arr.push(_s.value); if (i && _arr.length === i) break; } } catch (err) { _d = true; _e = err; } finally { try { if (!_n && _i["return"]) _i["return"](); } finally { if (_d) throw _e; } } return _arr; } return function (arr, i) { if (Array.isArray(arr)) { return arr; } else if (Symbol.iterator in Object(arr)) { return sliceIterator(arr, i); } else { throw new TypeError("Invalid attempt to destructure non-iterable instance"); } }; }(); /* eslint-disable no-labels */
/* eslint-disable no-fallthrough */
// SUBROUTINE ZASYI(ZR, ZI, FNU, KODE, N, YR, YI, NZ, RL, TOL, ELIM, ALIM)
// ***BEGIN PROLOGUE ZASYI
// ***REFER TO ZBESI,ZBESK
//
// ZASYI COMPUTES THE I BESSEL FUNCTION FOR REAL(Z).GE.0.0 BY
// MEANS OF THE ASYMPTOTIC EXPANSION FOR LARGE CABS(Z) IN THE
// REGION CABS(Z).GT.MAX(RL,FNU*FNU/2). NZ=0 IS A NORMAL RETURN.
// NZ.LT.0 INDICATES AN OVERFLOW ON KODE=1.
//
// ***ROUTINES CALLED D1MACH,AZABS,ZDIV,AZEXP,ZMLT,AZSQRT
// ***END PROLOGUE ZASYI


exports.zasyi = zasyi;

var _d1mach = require('../../utils/fortran-utils/d1mach.js');

var _zabs = require('./zabs.js');

var _zdiv3 = require('./zdiv.js');

var _zexp = require('./zexp.js');

var _zmlt7 = require('./zmlt.js');

var _zsqrt = require('./zsqrt.js');

function zasyi(zr, zi, fnu, kode, n, yr, yi, rl, tol, elim, alim) {
 var aa = void 0,
 aez = void 0,
 ak = void 0,
 ak1i = void 0,
 ak1r = void 0,
 arg = void 0,
 arm = void 0,
 atol = void 0,
 az = void 0,
 bb = void 0,
 bk = void 0,
 cki = void 0,
 ckr = void 0,
 conei = void 0,
 coner = void 0,
 cs1i = void 0,
 cs1r = void 0,
 cs2i = void 0,
 cs2r = void 0,
 czi = void 0,
 czr = void 0,
 dfnu = void 0,
 dki = void 0,
 dkr = void 0,
 dnu2 = void 0,
 ezi = void 0,
 ezr = void 0,
 fdn = void 0,
 pi = void 0,
 p1i = void 0,
 p1r = void 0,
 raz = void 0,
 rtpi = void 0,
 rtr1 = void 0,
 rzi = void 0,
 rzr = void 0,
 s = void 0,
 sgn = void 0,
 sqk = void 0,
 sti = void 0,
 str = void 0,
 s2i = void 0,
 s2r = void 0,
 tzi = void 0,
 tzr = void 0,
 zeroi = void 0,
 zeror = void 0,
 i = void 0,
 ib = void 0,
 il = void 0,
 inu = void 0,
 j = void 0,
 jl = void 0,
 k = void 0,
 koded = void 0,
 m = void 0,
 nn = void 0,
 nz = void 0;

 var goToLabel = 0;
 mainExecutionLoop: while (true) {
 mainSwitch: switch (goToLabel) {
 case 0:
 pi = 3.14159265358979324;
 rtpi = 0.159154943091895336;
 zeror = 0.0;
 zeroi = 0.0;
 coner = 1.0;
 conei = 0.0;

 nz = 0;
 az = (0, _zabs.azabs)(zr, zi);
 arm = 1.0e+3 * (0, _d1mach.d1mach)(1);
 rtr1 = Math.sqrt(arm);
 il = Math.min(2, n);
 dfnu = fnu + (n - il);
 // c-----------------------------------------------------------------------
 // c overflow test
 // c-----------------------------------------------------------------------
 raz = 1.0 / az;
 str = zr * raz;
 sti = -zi * raz;
 ak1r = rtpi * str * raz;
 ak1i = rtpi * sti * raz;

 var _azsqrt = (0, _zsqrt.azsqrt)(ak1r, ak1i);

 var _azsqrt2 = _slicedToArray(_azsqrt, 2);

 ak1r = _azsqrt2[0];
 ak1i = _azsqrt2[1];

 czr = zr;
 czi = zi;
 if (kode !== 2) {
 goToLabel = 10;break;
 }
 czr = zeror;
 czi = zi;
 case 10:
 if (Math.abs(czr) > elim) {
 goToLabel = 100;break;
 }
 dnu2 = dfnu + dfnu;
 koded = 1;
 if (Math.abs(czr) > alim && n > 2) {
 goToLabel = 20;break;
 }
 koded = 0;

 var _azexp = (0, _zexp.azexp)(czr, czi);

 var _azexp2 = _slicedToArray(_azexp, 2);

 str = _azexp2[0];
 sti = _azexp2[1];

 var _zmlt = (0, _zmlt7.zmlt)(ak1r, ak1i, str, sti);

 var _zmlt2 = _slicedToArray(_zmlt, 2);

 ak1r = _zmlt2[0];
 ak1i = _zmlt2[1];

 case 20:
 fdn = 0.0;
 if (dnu2 > rtr1) fdn = dnu2 * dnu2;
 ezr = zr * 8.0;
 ezi = zi * 8.0;
 // c-----------------------------------------------------------------------
 // c when z is imaginary, the error test must be made relative to the
 // c first reciprocal power since this is the leading term of the
 // c expansion for the imaginary part.
 // c-----------------------------------------------------------------------
 aez = 8.0 * az;
 s = tol / aez;
 jl = Math.trunc(rl + rl) + 2;
 p1r = zeror;
 p1i = zeroi;
 if (zi === 0.0) {
 goToLabel = 30;break;
 }
 // c-----------------------------------------------------------------------
 // c calculate exp(pi*(0.5+fnu+n-il)*i) to minimize losses of
 // c significance when fnu or n is large
 // c-----------------------------------------------------------------------
 inu = Math.trunc(fnu);
 arg = (fnu - inu) * pi;
 inu = inu + n - il;
 ak = -Math.sin(arg);
 bk = Math.cos(arg);
 if (zi < 0.0) {
 bk = -bk;
 }
 p1r = ak;
 p1i = bk;
 if (inu % 2 === 0) {
 goToLabel = 30;break;
 }
 p1r = -p1r;
 p1i = -p1i;
 case 30:
 // for loop 70:
 for (k = 1; k <= il; k++) {
 sqk = fdn - 1.0;
 atol = s * Math.abs(sqk);
 sgn = 1.0;
 cs1r = coner;
 cs1i = conei;
 cs2r = coner;
 cs2i = conei;
 ckr = coner;
 cki = conei;
 ak = 0.0;
 aa = 1.0;
 bb = aez;
 dkr = ezr;
 dki = ezi;
 // for loop 40
 for (j = 1; j <= jl; j++) {
 var _zdiv = (0, _zdiv3.zdiv)(ckr, cki, dkr, dki);

 var _zdiv2 = _slicedToArray(_zdiv, 2);

 str = _zdiv2[0];
 sti = _zdiv2[1];

 ckr = str * sqk;
 cki = sti * sqk;
 cs2r = cs2r + ckr;
 cs2i = cs2i + cki;
 sgn = -sgn;
 cs1r = cs1r + ckr * sgn;
 cs1i = cs1i + cki * sgn;
 dkr = dkr + ezr;
 dki = dki + ezi;
 aa = aa * Math.abs(sqk) / bb;
 bb = bb + aez;
 ak = ak + 8.0;
 sqk = sqk - ak;
 if (aa <= atol) {
 goToLabel = 50;break;
 }
 }
 if (goToLabel === 50) {
 // go to 50 - loop converged under atol
 } else {
 // 40 continue - throw error
 goToLabel = 110;break mainSwitch;
 }
 // 50 continue
 s2r = cs1r;
 s2i = cs1i;
 if (zr + zr >= elim) {
 // goToLabel = 60;
 } else {
 tzr = zr + zr;
 tzi = zi + zi;

 var _azexp3 = (0, _zexp.azexp)(-tzr, -tzi);

 var _azexp4 = _slicedToArray(_azexp3, 2);

 str = _azexp4[0];
 sti = _azexp4[1];

 var _zmlt3 = (0, _zmlt7.zmlt)(str, sti, p1r, p1i);

 var _zmlt4 = _slicedToArray(_zmlt3, 2);

 str = _zmlt4[0];
 sti = _zmlt4[1];

 var _zmlt5 = (0, _zmlt7.zmlt)(str, sti, cs2r, cs2i);

 var _zmlt6 = _slicedToArray(_zmlt5, 2);

 str = _zmlt6[0];
 sti = _zmlt6[1];

 s2r = s2r + str;
 s2i = s2i + sti;
 }
 fdn = fdn + 8.0 * dfnu + 4.0;
 p1r = -p1r;
 p1i = -p1i;
 m = n - il + k;
 yr[m - 1] = s2r * ak1r - s2i * ak1i;
 yi[m - 1] = s2r * ak1i + s2i * ak1r;
 }
 if (n <= 2) {
 break mainExecutionLoop;
 }
 nn = n;
 k = nn - 2;
 ak = k;
 str = zr * raz;
 sti = -zi * raz;
 rzr = (str + str) * raz;
 rzi = (sti + sti) * raz;
 ib = 3;
 // do 80 i=ib,nn
 for (i = ib; i <= nn; i++) {
 yr[k - 1] = (ak + fnu) * (rzr * yr[k] - rzi * yi[k]) + yr[k + 1];
 yi[k - 1] = (ak + fnu) * (rzr * yi[k] + rzi * yr[k]) + yi[k + 1];
 ak = ak - 1.0;
 k = k - 1;
 }
 if (koded === 0) {
 break mainExecutionLoop;
 }

 // do 90 i=1,nn
 var _azexp5 = (0, _zexp.azexp)(czr, czi);

 var _azexp6 = _slicedToArray(_azexp5, 2);

 ckr = _azexp6[0];
 cki = _azexp6[1];
 for (i = 1; i <= nn; i++) {
 str = yr[i - 1] * ckr - yi[i - 1] * cki;
 yi[i - 1] = yr[i - 1] * cki + yi[i - 1] * ckr;
 yr[i - 1] = str;
 }
 break mainExecutionLoop;
 case 100:
 nz = -1;
 break mainExecutionLoop;
 case 110:
 nz = -2;
 default:
 break mainExecutionLoop;
 }
 }

 return nz;
}
},{"../../utils/fortran-utils/d1mach.js":91,"./zabs.js":11,"./zdiv.js":26,"./zexp.js":27,"./zmlt.js":31,"./zsqrt.js":36}],16:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.zbesh = zbesh;

var _fortranHelpers = require('../../utils/fortranHelpers.js');

var fortranHelpers = _interopRequireWildcard(_fortranHelpers);

var _d1mach = require('../../utils/fortran-utils/d1mach.js');

var _i1mach = require('../../utils/fortran-utils/i1mach.js');

var _zuoik = require('./zuoik.js');

var _zbknu = require('./zbknu.js');

var _zacon = require('./zacon.js');

var _zbunk = require('./zbunk.js');

var _zabs = require('./zabs.js');

function _interopRequireWildcard(obj) { if (obj && obj.__esModule) { return obj; } else { var newObj = {}; if (obj != null) { for (var key in obj) { if (Object.prototype.hasOwnProperty.call(obj, key)) newObj[key] = obj[key]; } } newObj.default = obj; return newObj; } }

/* eslint-disable no-labels */
/* eslint-disable no-fallthrough */
// **BEGIN PROLOGUE ZBESH
// **DATE WRITTEN 830501 (YYMMDD)
// **REVISION DATE 890801 (YYMMDD)
// ***PORT TO ECMASCRIPT 201801 (YYYYMM)
// **CATEGORY NO. B5K
// **KEYWORDS H-BESSEL FUNCTIONS,BESSEL FUNCTIONS OF COMPLEX ARGUMENT,
// BESSEL FUNCTIONS OF THIRD KIND,HANKEL FUNCTIONS
// **AUTHOR (FORTRAN) AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
// **AUTHOR (ECMASCRIPT) ERB, KC, KINGS DISTRIBUTED SYSTEMS
// **PURPOSE TO COMPUTE THE H-BESSEL FUNCTIONS OF A COMPLEX ARGUMENT
// **DESCRIPTION
//
// ***A DOUBLE PRECISION ROUTINE***
// ON KODE=1, ZBESH COMPUTES AN N MEMBER SEQUENCE OF COMPLEX
// HANKEL (BESSEL) FUNCTIONS CY(J)=H(M,FNU+J-1,Z) FOR KINDS M=1
// OR 2, REAL, NONNEGATIVE ORDERS FNU+J-1, J=1,...,N, AND COMPLEX
// Z !== CMPLX(0.0,0.0) IN THE CUT PLANE -PI < ARG(Z) <= PI.
// ON KODE=2, ZBESH RETURNS THE SCALED HANKEL FUNCTIONS
//
// CY(I)=MATH.EXP(-MM*Z*I)*H(M,FNU+J-1,Z) MM=3-2*M, I**2=-1.
//
// WHICH REMOVES THE EXPONENTIAL BEHAVIOR IN BOTH THE UPPER AND
// LOWER HALF PLANES. DEFINITIONS AND NOTATION ARE FOUND IN THE
// NBS HANDBOOK OF MATHEMATICAL FUNCTIONS (REF. 1).
//
// INPUT ZR,ZI,FNU ARE DOUBLE PRECISION
// ZR,ZI - Z=CMPLX(ZR,ZI), Z !== CMPLX(0.0E0,0.0E0),
// -PT < ARG(Z) <= PI
// FNU - ORDER OF INITIAL H FUNCTION, FNU >= 0.0E0
// KODE - A PARAMETER TO INDICATE THE SCALING OPTION
// KODE= 1 RETURNS
// CY(J)=H(M,FNU+J-1,Z), J=1,...,N
// = 2 RETURNS
// CY(J)=H(M,FNU+J-1,Z)*MATH.EXP(-I*Z*(3-2M))
// J=1,...,N , I**2=-1
// M - KIND OF HANKEL FUNCTION, M=1 OR 2
// N - NUMBER OF MEMBERS IN THE SEQUENCE, N >= 1
//
// OUTPUT CYR,CYI ARE DOUBLE PRECISION
// CYR,CYI- DOUBLE PRECISION VECTORS WHOSE FIRST N COMPONENTS
// CONTAIN REAL AND IMAGINARY PARTS FOR THE SEQUENCE
// CY(J)=H(M,FNU+J-1,Z) OR
// CY(J)=H(M,FNU+J-1,Z)*MATH.EXP(-I*Z*(3-2M)) J=1,...,N
// DEPENDING ON KODE, I**2=-1.
// NZ - NUMBER OF COMPONENTS SET TO ZERO DUE TO UNDERFLOW,
// NZ= 0 , NORMAL RETURN
// NZ > 0 , FIRST NZ COMPONENTS OF CY SET TO ZERO DUE
// TO UNDERFLOW, CY(J)=CMPLX(0.0E0,0.0E0)
// J=1,...,NZ WHEN Y > 0.0 AND M=1 OR
// Y < 0.0 AND M=2. FOR THE COMPLMENTARY
// HALF PLANES, NZ STATES ONLY THE NUMBER
// OF UNDERFLOWS.
// IERR - ERROR FLAG
// IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
// IERR=1, INPUT ERROR - NO COMPUTATION
// IERR=2, OVERFLOW - NO COMPUTATION, FNU TOO
// LARGE OR CABS(Z) TOO SMALL OR BOTH
// IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE
// BUT LOSSES OF SIGNIFCANCE BY ARGUMENT
// REDUCTION PRODUCE LESS THAN HALF OF MACHINE
// ACCURACY
// IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTA-
// TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-
// CANCE BY ARGUMENT REDUCTION
// IERR=5, ERROR - NO COMPUTATION,
// ALGORITHM TERMINATION CONDITION NOT MET
//
// **LONG DESCRIPTION
//
// THE COMPUTATION IS CARRIED OUT BY THE RELATION
//
// H(M,FNU,Z)=(MATH.TRUNC(1/M)P)*MATH.EXP(-MP*FNU)*K(FNU,Z*MATH.EXP(-MP))
// MP=MM*HPI*I, MM=3-2*M, HPI=PMATH.TRUNC(I/2), I**2=-1
//
// FOR M=1 OR 2 WHERE THE K BESSEL FUNCTION IS COMPUTED FOR THE
// RIGHT HALF PLANE RE(Z) >= 0.0. THE K FUNCTION IS CONTINUED
// TO THE LEFT HALF PLANE BY THE RELATION
//
// K(FNU,Z*MATH.EXP(MP)) = MATH.EXP(-MP*FNU)*K(FNU,Z)-MP*I(FNU,Z)
// MP=MR*PI*I, MR=+1 OR -1, RE(Z) > 0, I**2=-1
//
// WHERE I(FNU,Z) IS THE I BESSEL FUNCTION.
//
// EXPONENTIAL DECAY OF H(M,FNU,Z) OCCURS IN THE UPPER HALF Z
// PLANE FOR M=1 AND THE LOWER HALF Z PLANE FOR M=2. EXPONENTIAL
// GROWTH OCCURS IN THE COMPLEMENTARY HALF PLANES. SCALING
// BY MATH.EXP(-MM*Z*I) REMOVES THE EXPONENTIAL BEHAVIOR IN THE
// WHOLE Z PLANE FOR Z TO INFINITY.
//
// FOR NEGATIVE ORDERS,THE FORMULAE
//
// H(1,-FNU,Z) = H(1,FNU,Z)*CEXP( PI*FNU*I)
// H(2,-FNU,Z) = H(2,FNU,Z)*CEXP(-PI*FNU*I)
// I**2=-1
//
// CAN BE USED.
//
// IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
// MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z OR FNU+N-1 IS
// LARGE, LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR.
// CONSEQUENTLY, IF EITHER ONE EXCEEDS U1=MATH.SQRT(0.5/UR), {
// LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR FLAG
// IERR=3 IS TRIGGERED WHERE UR=MATH.MAX(D1MACH(4),1.0E-18) IS
// DOUBLE PRECISION UNIT ROUNDOFF LIMITED TO 18 DIGITS PRECISION.
// IF EITHER IS LARGER THAN U2=0.5/UR, { ALL SIGNIFICANCE IS
// LOST AND IERR=4. IN ORDER TO USE THE INT FUNCTION, ARGUMENTS
// MUST BE FURTHER RESTRICTED NOT TO EXCEED THE LARGEST MACHINE
// INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF Z AND FNU+N-1 IS
// RESTRICTED BY MATH.MIN(U2,U3). ON 32 BIT MACHINES, U1,U2, AND U3
// ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE PRECISION
// ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE PRECISION
// ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMITING IN
// THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT ONE CAN EXPECT
// TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES, NO DIGITS
// IN SINGLE AND ONLY 7 DIGITS IN DOUBLE PRECISION ARITHMETIC.
// SIMILAR CONSIDERATIONS HOLD FOR OTHER MACHINES.
//
// THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
// BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MATH.MAX(UNIT
// ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
// SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
// ELEMENTARY FUNCTIONS. HERE, S=MATH.MAX(1,MATH.ABS(MATH.LOG10(CABS(Z))),
// MATH.ABS(MATH.LOG10(FNU))) APPROXIMATELY (I.E. S=MATH.MAX(1,MATH.ABS(EXPONENT OF
// CABS(Z),ABS(EXPONENT OF FNU))) ). HOWEVER, THE PHASE ANGLE MAY
// HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
// ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
// SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
// THAN THE OTHER, { ONE CAN EXPECT ONLY MATH.MAX(MATH.ABS(MATH.LOG10(P))-K,
// 0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
// THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
// COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
// BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
// COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
// MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
// THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PMATH.TRUNC(I/2)-P,
// OR -PMATH.TRUNC(I/2)+P.
//
// **REFERENCES HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
// AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
// COMMERCE, 1955.
//
// COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
// BY D. E. AMOS, SAND83-0083, MAY, 1983.
//
// COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
// AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
//
// A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
// ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
// 1018, MAY, 1985
//
// A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
// ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, TRANS.
// MATH. SOFTWARE, 1986
//
// **ROUTINES ED() ZACON,ZBKNU,ZBUNK,ZUOIK,AZABS,I1MACH,D1MACH
// **END PROLOGUE ZBESH
//
// COMPLEX CY,Z,ZN,ZT,CSGN
function zbesh(zr, zi, fnu, kode, m, n, cyr, cyi) {
 var aa = void 0,
 alim = void 0,
 aln = void 0,
 arg = void 0,
 az = void 0,
 dig = void 0,
 elim = void 0,
 fmm = void 0,
 fn = void 0,
 fnul = void 0,
 hpi = void 0,
 rhpi = void 0,
 rl = void 0,
 r1m5 = void 0,
 sgn = void 0,
 str = void 0,
 tol = void 0,
 ufl = void 0,
 zni = void 0,
 znr = void 0,
 zti = void 0,
 bb = void 0,
 ascle = void 0,
 rtol = void 0,
 atol = void 0,
 sti = void 0,
 csgnr = void 0,
 csgni = void 0,
 i = void 0,
 inu = void 0,
 inuh = void 0,
 ir = void 0,
 k = void 0,
 k1 = void 0,
 k2 = void 0,
 mm = void 0,
 mr = void 0,
 nn = void 0,
 nuf = void 0,
 nw = void 0,
 nz = void 0,
 ierr = void 0;

 var goToLabel = 0;
 mainExecutionLoop: while (true) {
 switch (goToLabel) {
 case 0:
 hpi = 1.57079632679489662e0;
 ierr = 0;
 nz = 0;
 if (zr === 0.0e0 && zi === 0.0e0) ierr = 1;
 if (fnu < 0.0e0) ierr = 1;
 if (m < 1 || m > 2) ierr = 1;
 if (kode < 1 || kode > 2) ierr = 1;
 if (n < 1) ierr = 1;
 if (ierr !== 0) break mainExecutionLoop;
 nn = n;
 // -----------------------------------------------------------------------
 // set parameters related to machine constants.
 // tol is the approximate unit roundoff limited to 1.0e-18.
 // elim is the approximate exponential over- and underflow limit.
 // Math.exp(-elim) < Math.exp(-alim)=Math.exp(-elim)/tol and
 // Math.exp(elim) > Math.exp(alim)=Math.exp(elim)*tol are intervals near
 // underflow and overflow limits where scaled arithmetic is done.
 // rl is the lower boundary of the asymptotic expansion for large z.
 // dig = number of base 10 digits in tol = 10**(-dig).
 // fnul is the lower boundary of the asymptotic series for large fnu
 // -----------------------------------------------------------------------
 tol = Math.max((0, _d1mach.d1mach)(4), 1.0e-18);
 k1 = (0, _i1mach.i1mach)(15);
 k2 = (0, _i1mach.i1mach)(16);
 r1m5 = (0, _d1mach.d1mach)(5);
 k = Math.min(Math.abs(k1), Math.abs(k2));
 elim = 2.303e0 * (k * r1m5 - 3.0e0);
 k1 = (0, _i1mach.i1mach)(14) - 1;
 aa = r1m5 * k1;
 dig = Math.min(aa, 18.0e0);
 aa = aa * 2.303e0;
 alim = elim + Math.max(-aa, -41.45e0);
 fnul = 10.0e0 + 6.0e0 * (dig - 3.0e0);
 rl = 1.2e0 * dig + 3.0e0;
 fn = fnu + (nn - 1);
 mm = 3 - m - m;
 fmm = mm;
 znr = fmm * zi;
 zni = -fmm * zr;
 // -----------------------------------------------------------------------
 // test for proper range
 // -----------------------------------------------------------------------
 az = (0, _zabs.azabs)(zr, zi);
 aa = 0.5e0 / tol;
 bb = (0, _i1mach.i1mach)(9) * 0.5e0;
 aa = Math.min(aa, bb);
 if (az > aa) {
 goToLabel = 260;break;
 }
 if (fn > aa) {
 goToLabel = 260;break;
 }
 aa = Math.sqrt(aa);
 if (az > aa) ierr = 3;
 if (fn > aa) ierr = 3;
 // -----------------------------------------------------------------------
 // overflow test on the last member of the sequence
 // -----------------------------------------------------------------------
 ufl = (0, _d1mach.d1mach)(1) * 1.0e+3;
 if (az < ufl) {
 goToLabel = 230;break;
 }
 if (fnu > fnul) {
 goToLabel = 90;break;
 }
 if (fn <= 1.0e0) {
 goToLabel = 70;break;
 }
 if (fn > 2.0e0) {
 goToLabel = 60;break;
 }
 if (az > tol) {
 goToLabel = 70;break;
 }
 arg = 0.5e0 * az;
 aln = -fn * Math.log(arg);
 if (aln > elim) {
 goToLabel = 230;break;
 }
 goToLabel = 70;break;
 case 60:

 nuf = (0, _zuoik.zuoik)(znr, zni, fnu, kode, 2, nn, cyr, cyi, tol, elim, alim);
 if (nuf < 0) {
 goToLabel = 230;break;
 }
 nz = nz + nuf;
 nn = nn - nuf;
 // -----------------------------------------------------------------------
 // here nn=n or nn=0 since nuf=0,nn, or -1 on return from cuoik
 // if nuf=nn, { cy(i)=czero for all i
 // -----------------------------------------------------------------------
 if (nn === 0) {
 goToLabel = 140;break;
 }
 case 70:

 if (znr < 0.0e0 || znr === 0.0e0 && zni < 0.0e0 && m === 2) {
 goToLabel = 80;break;
 }
 // -----------------------------------------------------------------------
 // right half plane computation, xn >= 0. && (xn !== 0. ||
 // yn >= 0. || m=1)
 // -----------------------------------------------------------------------
 nz = (0, _zbknu.zbknu)(znr, zni, fnu, kode, nn, cyr, cyi, tol, elim, alim);
 goToLabel = 110;break;
 // -----------------------------------------------------------------------
 // left half plane computation
 // -----------------------------------------------------------------------
 case 80:

 mr = -mm;
 nw = (0, _zacon.zacon)(znr, zni, fnu, kode, mr, nn, cyr, cyi, rl, fnul, tol, elim, alim);
 if (nw < 0) {
 goToLabel = 240;break;
 }
 nz = nw;
 goToLabel = 110;break;
 case 90:

 // -----------------------------------------------------------------------
 // uniform asymptotic expansions for fnu > fnul
 // -----------------------------------------------------------------------
 mr = 0;
 if (znr >= 0.0e0 && (znr !== 0.0e0 || zni >= 0.0e0 || m !== 2)) {
 goToLabel = 100;break;
 }
 mr = -mm;
 if (znr !== 0.0e0 || zni >= 0.0e0) {
 goToLabel = 100;break;
 }
 znr = -znr;
 zni = -zni;
 case 100:

 nw = (0, _zbunk.zbunk)(znr, zni, fnu, kode, mr, nn, cyr, cyi, nw, tol, elim, alim);
 if (nw < 0) {
 goToLabel = 240;break;
 }
 nz = nz + nw;
 case 110:

 // -----------------------------------------------------------------------
 // h(m,fnu,z) = -fmm*(i/hpi)*(zt**fnu)*k(fnu,-z*zt)
 //
 // zt=Math.exp(-fmm*hpi*i) = cmplx(0.0,-fmm), fmm=3-2*m, m=1,2
 // -----------------------------------------------------------------------
 sgn = fortranHelpers.sign(hpi, -fmm);
 // -----------------------------------------------------------------------
 // calculate Math.exp(fnu*hpi*i) to minimize losses of significance
 // when fnu is large
 // -----------------------------------------------------------------------
 inu = Math.trunc(fnu);
 inuh = Math.trunc(inu / 2);
 ir = inu - 2 * inuh;
 arg = (fnu - (inu - ir)) * sgn;
 rhpi = 1.0e0 / sgn;
 // zni = rhpi*Math.cos(arg)
 // znr = -rhpi*Math.sin(arg)
 csgni = rhpi * Math.cos(arg);
 csgnr = -rhpi * Math.sin(arg);
 if (inuh % 2 === 0) {
 goToLabel = 120;break;
 }
 // znr = -znr
 // zni = -zni
 csgnr = -csgnr;
 csgni = -csgni;
 case 120:

 zti = -fmm;
 rtol = 1.0e0 / tol;
 ascle = ufl * rtol;
 for (i = 1; i <= nn; i++) {
 // str = cyr(i)*znr - cyi(i)*zni
 // cyi(i) = cyr(i)*zni + cyi(i)*znr
 // cyr(i) = str
 // str = -zni*zti
 // zni = znr*zti
 // znr = str
 aa = cyr[i - 1];
 bb = cyi[i - 1];
 atol = 1.0e0;
 if (Math.max(Math.abs(aa), Math.abs(bb)) > ascle) {
 // goToLabel = 135; break;
 } else {
 aa = aa * rtol;
 bb = bb * rtol;
 atol = tol;
 }
 // case 135:
 str = aa * csgnr - bb * csgni;
 sti = aa * csgni + bb * csgnr;
 cyr[i - 1] = str * atol;
 cyi[i - 1] = sti * atol;
 str = -csgni * zti;
 csgni = csgnr * zti;
 csgnr = str;
 }
 break mainExecutionLoop;
 case 140:

 if (znr < 0.0e0) {
 goToLabel = 230;break;
 }
 break mainExecutionLoop;
 case 230:

 nz = 0;
 ierr = 2;
 break mainExecutionLoop;
 case 240:

 if (nw === -1) {
 goToLabel = 230;break;
 }
 nz = 0;
 ierr = 5;
 break mainExecutionLoop;
 case 260:

 nz = 0;
 ierr = 4;

 default:
 break mainExecutionLoop;
 }
 }
 return [nz, ierr];
}
},{"../../utils/fortran-utils/d1mach.js":91,"../../utils/fortran-utils/i1mach.js":92,"../../utils/fortranHelpers.js":93,"./zabs.js":11,"./zacon.js":13,"./zbknu.js":23,"./zbunk.js":25,"./zuoik.js":44}],17:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.zbesi = zbesi;

var _i1mach = require('../../utils/fortran-utils/i1mach.js');

var _d1mach = require('../../utils/fortran-utils/d1mach.js');

var _zbinu = require('./zbinu.js');

var _zabs = require('./zabs.js');

/* eslint-disable no-labels */
/* eslint-disable no-fallthrough */
// **BEGIN PROLOGUE ZBESI
// **DATE WRITTEN 830501 (YYMMDD)
// **REVISION DATE 890801 (YYMMDD)
// ***PORT TO ECMASCRIPT 201801 (YYYYMM)
// **CATEGORY NO. B5K
// **KEYWORDS I-BESSEL FUNCTION,COMPLEX BESSEL FUNCTION,
// MODIFIED BESSEL FUNCTION OF THE FIRST KIND
// **AUTHOR (FORTRAN) AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
// **AUTHOR (ECMASCRIPT) ERB, KC, KINGS DISTRIBUTED SYSTEMS
// **PURPOSE TO COMPUTE I-BESSEL FUNCTIONS OF COMPLEX ARGUMENT
// **DESCRIPTION
//
// ***A DOUBLE PRECISION ROUTINE***
// ON KODE=1, ZBESI COMPUTES AN N MEMBER SEQUENCE OF COMPLEX
// BESSEL FUNCTIONS CY(J)=I(FNU+J-1,Z) FOR REAL, NONNEGATIVE
// ORDERS FNU+J-1, J=1,...,N AND COMPLEX Z IN THE CUT PLANE
// -PI < ARG(Z) <= PI. ON KODE=2, ZBESI RETURNS THE SCALED
// FUNCTIONS
//
// CY(J)=MATH.EXP(-MATH.ABS(X))*I(FNU+J-1,Z) J = 1,...,N , X=COMPLEXHELPERS.RE(Z)
//
// WITH THE EXPONENTIAL GROWTH REMOVED IN BOTH THE LEFT AND
// RIGHT HALF PLANES FOR Z TO INFINITY. DEFINITIONS AND NOTATION
// ARE FOUND IN THE NBS HANDBOOK OF MATHEMATICAL FUNCTIONS
// (REF. 1).
//
// INPUT ZR,ZI,FNU ARE DOUBLE PRECISION
// ZR,ZI - Z=CMPLX(ZR,ZI), -PI < ARG(Z) <= PI
// FNU - ORDER OF INITIAL I FUNCTION, FNU >= 0.0E0
// KODE - A PARAMETER TO INDICATE THE SCALING OPTION
// KODE= 1 RETURNS
// CY(J)=I(FNU+J-1,Z), J=1,...,N
// = 2 RETURNS
// CY(J)=I(FNU+J-1,Z)*MATH.EXP(-MATH.ABS(X)), J=1,...,N
// N - NUMBER OF MEMBERS OF THE SEQUENCE, N >= 1
//
// OUTPUT CYR,CYI ARE DOUBLE PRECISION
// CYR,CYI- DOUBLE PRECISION VECTORS WHOSE FIRST N COMPONENTS
// CONTAIN REAL AND IMAGINARY PARTS FOR THE SEQUENCE
// CY(J)=I(FNU+J-1,Z) OR
// CY(J)=I(FNU+J-1,Z)*MATH.EXP(-MATH.ABS(X)) J=1,...,N
// DEPENDING ON KODE, X=COMPLEXHELPERS.RE(Z)
// NZ - NUMBER OF COMPONENTS SET TO ZERO DUE TO UNDERFLOW,
// NZ= 0 , NORMAL RETURN
// NZ > 0 , LAST NZ COMPONENTS OF CY SET TO ZERO
// TO UNDERFLOW, CY(J)=CMPLX(0.0E0,0.0E0)
// J = N-NZ+1,...,N
// IERR - ERROR FLAG
// IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
// IERR=1, INPUT ERROR - NO COMPUTATION
// IERR=2, OVERFLOW - NO COMPUTATION, RE(Z) TOO
// LARGE ON KODE=1
// IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE
// BUT LOSSES OF SIGNIFCANCE BY ARGUMENT
// REDUCTION PRODUCE LESS THAN HALF OF MACHINE
// ACCURACY
// IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTA-
// TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-
// CANCE BY ARGUMENT REDUCTION
// IERR=5, ERROR - NO COMPUTATION,
// ALGORITHM TERMINATION CONDITION NOT MET
//
// **LONG DESCRIPTION
//
// THE COMPUTATION IS CARRIED OUT BY THE POWER SERIES FOR
// SMALL CABS(Z), THE ASYMPTOTIC EXPANSION FOR LARGE CABS(Z),
// THE MILLER ALGORITHM NORMALIZED BY THE WRONSKIAN AND A
// NEUMANN SERIES FOR IMTERMEDIATE MAGNITUDES, AND THE
// UNIFORM ASYMPTOTIC EXPANSIONS FOR I(FNU,Z) AND J(FNU,Z)
// FOR LARGE ORDERS. BACKWARD RECURRENCE IS USED TO GENERATE
// SEQUENCES OR REDUCE ORDERS WHEN NECESSARY.
//
// THE CALCULATIONS ABOVE ARE DONE IN THE RIGHT HALF PLANE AND
// CONTINUED INTO THE LEFT HALF PLANE BY THE FORMULA
//
// I(FNU,Z*MATH.EXP(M*PI)) = MATH.EXP(M*PI*FNU)*I(FNU,Z) COMPLEXHELPERS.RE(Z) > 0.0
// M = +I OR -I, I**2=-1
//
// FOR NEGATIVE ORDERS,THE FORMULA
//
// I(-FNU,Z) = I(FNU,Z) + (2/PI)*MATH.SIN(PI*FNU)*K(FNU,Z)
//
// CAN BE USED. HOWEVER,FOR LARGE ORDERS CLOSE TO INTEGERS, THE
// THE FUNCTION CHANGES RADIY(). WHEN FNU IS A LARGE POSITIVE
// INTEGER,THE MAGNITUDE OF I(-FNU,Z)=I(FNU,Z) IS A LARGE
// NEGATIVE POWER OF TEN. BUT WHEN FNU IS NOT AN INTEGER,
// K(FNU,Z) DOMINATES IN MAGNITUDE WITH A LARGE POSITIVE POWER OF
// TEN AND THE MOST THAT THE SECOND TERM CAN BE REDUCED IS BY
// UNIT ROUNDOFF FROM THE COEFFICIENT. THUS, WIDE CHANGES CAN
// OCCUR WITHIN UNIT ROUNDOFF OF A LARGE INTEGER FOR FNU. HERE,
// LARGE MEANS FNU > CABS(Z).
//
// IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
// MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z OR FNU+N-1 IS
// LARGE, LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR.
// CONSEQUENTLY, IF EITHER ONE EXCEEDS U1=MATH.SQRT(0.5/UR), {
// LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR FLAG
// IERR=3 IS TRIGGERED WHERE UR=MATH.MAX(D1MACH(4),1.0E-18) IS
// DOUBLE PRECISION UNIT ROUNDOFF LIMITED TO 18 DIGITS PRECISION.
// IF EITHER IS LARGER THAN U2=0.5/UR, { ALL SIGNIFICANCE IS
// LOST AND IERR=4. IN ORDER TO USE THE INT FUNCTION, ARGUMENTS
// MUST BE FURTHER RESTRICTED NOT TO EXCEED THE LARGEST MACHINE
// INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF Z AND FNU+N-1 IS
// RESTRICTED BY MATH.MIN(U2,U3). ON 32 BIT MACHINES, U1,U2, AND U3
// ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE PRECISION
// ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE PRECISION
// ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMITING IN
// THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT ONE CAN EXPECT
// TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES, NO DIGITS
// IN SINGLE AND ONLY 7 DIGITS IN DOUBLE PRECISION ARITHMETIC.
// SIMILAR CONSIDERATIONS HOLD FOR OTHER MACHINES.
//
// THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
// BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MATH.MAX(UNIT
// ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
// SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
// ELEMENTARY FUNCTIONS. HERE, S=MATH.MAX(1,MATH.ABS(MATH.LOG10(CABS(Z))),
// MATH.ABS(MATH.LOG10(FNU))) APPROXIMATELY (I.E. S=MATH.MAX(1,MATH.ABS(EXPONENT OF
// CABS(Z),ABS(EXPONENT OF FNU))) ). HOWEVER, THE PHASE ANGLE MAY
// HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
// ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
// SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
// THAN THE OTHER, { ONE CAN EXPECT ONLY MATH.MAX(MATH.ABS(MATH.LOG10(P))-K,
// 0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
// THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
// COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
// BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
// COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
// MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
// THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PMATH.TRUNC(I/2)-P,
// OR -PMATH.TRUNC(I/2)+P.
//
// **REFERENCES HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
// AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
// COMMERCE, 1955.
//
// COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
// BY D. E. AMOS, SAND83-0083, MAY, 1983.
//
// COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
// AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
//
// A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
// ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
// 1018, MAY, 1985
//
// A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
// ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, TRANS.
// MATH. SOFTWARE, 1986
//
// **ROUTINES ED() ZBINU,I1MACH,D1MACH
// **END PROLOGUE ZBESI
// COMPLEX CONE,CSGN,CW,CY,CZERO,Z,ZN
function zbesi(zr, zi, fnu, kode, n, cyr, cyi) {
 var aa = void 0,
 alim = void 0,
 arg = void 0,
 conei = void 0,
 coner = void 0,
 csgni = void 0,
 csgnr = void 0,
 dig = void 0,
 elim = void 0,
 fnul = void 0,
 pi = void 0,
 rl = void 0,
 r1m5 = void 0,
 str = void 0,
 tol = void 0,
 zni = void 0,
 znr = void 0,
 az = void 0,
 bb = void 0,
 fn = void 0,
 ascle = void 0,
 rtol = void 0,
 atol = void 0,
 sti = void 0,
 i = void 0,
 inu = void 0,
 k = void 0,
 k1 = void 0,
 k2 = void 0,
 nn = void 0,
 ierr = void 0,
 nz = void 0;

 var goToLabel = 0;
 mainExecutionLoop: while (true) {
 switch (goToLabel) {
 case 0:
 pi = 3.14159265358979324e0;
 coner = 1.0e0;
 conei = 0.0e0;

 ierr = 0;
 nz = 0;
 if (fnu < 0.0e0) ierr = 1;
 if (kode < 1 || kode > 2) ierr = 1;
 if (n < 1) ierr = 1;
 if (ierr !== 0) break mainExecutionLoop;
 // -----------------------------------------------------------------------
 // set parameters related to machine constants.
 // tol is the approximate unit roundoff limited to 1.0e-18.
 // elim is the approximate exponential over- and underflow limit.
 // Math.exp(-elim) < Math.exp(-alim)=Math.exp(-elim)/tol and
 // Math.exp(elim) > Math.exp(alim)=Math.exp(elim)*tol are intervals near
 // underflow and overflow limits where scaled arithmetic is done.
 // rl is the lower boundary of the asymptotic expansion for large z.
 // dig = number of base 10 digits in tol = 10**(-dig).
 // fnul is the lower boundary of the asymptotic series for large fnu.
 // -----------------------------------------------------------------------
 tol = Math.max((0, _d1mach.d1mach)(4), 1.0e-18);
 k1 = (0, _i1mach.i1mach)(15);
 k2 = (0, _i1mach.i1mach)(16);
 r1m5 = (0, _d1mach.d1mach)(5);
 k = Math.min(Math.abs(k1), Math.abs(k2));
 elim = 2.303e0 * (k * r1m5 - 3.0e0);
 k1 = (0, _i1mach.i1mach)(14) - 1;
 aa = r1m5 * k1;
 dig = Math.min(aa, 18.0e0);
 aa = aa * 2.303e0;
 alim = elim + Math.max(-aa, -41.45e0);
 rl = 1.2e0 * dig + 3.0e0;
 fnul = 10.0e0 + 6.0e0 * (dig - 3.0e0);
 // -----------------------------------------------------------------------------
 // test for proper range
 // -----------------------------------------------------------------------
 az = (0, _zabs.azabs)(zr, zi);
 fn = fnu + (n - 1);
 aa = 0.5e0 / tol;
 bb = (0, _i1mach.i1mach)(9) * 0.5e0;
 aa = Math.min(aa, bb);
 if (az > aa) {
 goToLabel = 260;break;
 }
 if (fn > aa) {
 goToLabel = 260;break;
 }
 aa = Math.sqrt(aa);
 if (az > aa) ierr = 3;
 if (fn > aa) ierr = 3;
 znr = zr;
 zni = zi;
 csgnr = coner;
 csgni = conei;
 if (zr >= 0.0e0) {
 goToLabel = 40;break;
 }
 znr = -zr;
 zni = -zi;
 // -----------------------------------------------------------------------
 // calculate csgn=Math.exp(fnu*pi*i) to minimize losses of significance
 // when fnu is large
 // -----------------------------------------------------------------------
 inu = Math.trunc(fnu);
 arg = (fnu - inu) * pi;
 if (zi < 0.0e0) arg = -arg;
 csgnr = Math.cos(arg);
 csgni = Math.sin(arg);
 if (inu % 2 === 0) {
 goToLabel = 40;break;
 }
 csgnr = -csgnr;
 csgni = -csgni;
 case 40:

 nz = (0, _zbinu.zbinu)(znr, zni, fnu, kode, n, cyr, cyi, rl, fnul, tol, elim, alim);
 if (nz < 0) {
 goToLabel = 120;break;
 }
 if (zr >= 0.0e0) break mainExecutionLoop;
 // -----------------------------------------------------------------------
 // analytic continuation to the left half plane
 // -----------------------------------------------------------------------
 nn = n - nz;
 if (nn === 0) break mainExecutionLoop;
 rtol = 1.0e0 / tol;
 ascle = (0, _d1mach.d1mach)(1) * rtol * 1.0e+3;
 for (i = 1; i <= nn; i++) {
 // str = cyr(i)*csgnr - cyi(i)*csgni
 // cyi(i) = cyr(i)*csgni + cyi(i)*csgnr
 // cyr(i) = str
 aa = cyr[i - 1];
 bb = cyi[i - 1];
 atol = 1.0e0;
 if (Math.max(Math.abs(aa), Math.abs(bb)) > ascle) {
 // goToLabel = 55; break;
 } else {
 aa = aa * rtol;
 bb = bb * rtol;
 atol = tol;
 }
 // case 55:
 str = aa * csgnr - bb * csgni;
 sti = aa * csgni + bb * csgnr;
 cyr[i - 1] = str * atol;
 cyi[i - 1] = sti * atol;
 csgnr = -csgnr;
 csgni = -csgni;
 }
 break mainExecutionLoop;
 case 120:

 if (nz === -2) {
 goToLabel = 130;break;
 }
 nz = 0;
 ierr = 2;
 break mainExecutionLoop;
 case 130:

 nz = 0;
 ierr = 5;
 break mainExecutionLoop;
 case 260:

 nz = 0;
 ierr = 4;

 default:
 break mainExecutionLoop;
 }
 }
 return [nz, ierr];
}
},{"../../utils/fortran-utils/d1mach.js":91,"../../utils/fortran-utils/i1mach.js":92,"./zabs.js":11,"./zbinu.js":21}],18:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.zbesj = zbesj;

var _d1mach = require('../../utils/fortran-utils/d1mach.js');

var _i1mach = require('../../utils/fortran-utils/i1mach.js');

var _zbinu = require('./zbinu.js');

function zbesj(zr, zi, fnu, kode, n, cyr, cyi) {
 var aa = void 0,
 alim = void 0,
 arg = void 0,
 cii = void 0,
 csgni = void 0,
 csgnr = void 0,
 dig = void 0,
 elim = void 0,
 fnul = void 0,
 hpi = void 0,
 rl = void 0,
 r1m5 = void 0,
 str = void 0,
 tol = void 0,
 zni = void 0,
 znr = void 0,
 bb = void 0,
 fn = void 0,
 az = void 0,
 ascle = void 0,
 rtol = void 0,
 atol = void 0,
 sti = void 0,
 i = void 0,
 inu = void 0,
 inuh = void 0,
 ir = void 0,
 k = void 0,
 k1 = void 0,
 k2 = void 0,
 nl = void 0,
 nz = void 0,
 ierr = void 0;

 var goToLabel = 0;
 mainExecutionLoop: while (true) {
 switch (goToLabel) {
 case 0:
 hpi = 1.57079632679489662e0;
 ierr = 0;
 nz = 0;
 if (fnu < 0.0e0) ierr = 1;
 if (kode < 1 || kode > 2) ierr = 1;
 if (n < 1) ierr = 1;
 if (ierr !== 0) break mainExecutionLoop;
 // -----------------------------------------------------------------------
 // set parameters related to machine constants.
 // tol is the approximate unit roundoff limited to 1.0e-18.
 // elim is the approximate exponential over- and underflow limit.
 // Math.exp(-elim) < Math.exp(-alim)=Math.exp(-elim)/tol and
 // Math.exp(elim) > Math.exp(alim)=Math.exp(elim)*tol are intervals near
 // underflow and overflow limits where scaled arithmetic is done.
 // rl is the lower boundary of the asymptotic expansion for large z.
 // dig = number of base 10 digits in tol = 10**(-dig).
 // fnul is the lower boundary of the asymptotic series for large fnu.
 // -----------------------------------------------------------------------
 tol = Math.max((0, _d1mach.d1mach)(4), 1.0e-18);
 k1 = (0, _i1mach.i1mach)(15);
 k2 = (0, _i1mach.i1mach)(16);
 r1m5 = (0, _d1mach.d1mach)(5);
 k = Math.min(Math.abs(k1), Math.abs(k2));
 elim = 2.303e0 * (k * r1m5 - 3.0e0);
 k1 = (0, _i1mach.i1mach)(14) - 1;
 aa = r1m5 * k1;
 dig = Math.min(aa, 18.0e0);
 aa = aa * 2.303e0;
 alim = elim + Math.max(-aa, -41.45e0);
 rl = 1.2e0 * dig + 3.0e0;
 fnul = 10.0e0 + 6.0e0 * (dig - 3.0e0);
 // -----------------------------------------------------------------------
 // test for proper range
 // -----------------------------------------------------------------------
 az = Math.abs(zr, zi);
 fn = fnu + (n - 1);
 aa = 0.5e0 / tol;
 bb = (0, _i1mach.i1mach)(9) * 0.5e0;
 aa = Math.min(aa, bb);
 if (az > aa) {
 goToLabel = 260;break;
 }
 if (fn > aa) {
 goToLabel = 260;break;
 }
 aa = Math.sqrt(aa);
 if (az > aa) ierr = 3;
 if (fn > aa) ierr = 3;
 // -----------------------------------------------------------------------
 // calculate csgn=Math.exp(fnu*hpi*i) to minimize losses of significance
 // when fnu is large
 // -----------------------------------------------------------------------
 cii = 1.0e0;
 inu = Math.trunc(fnu);
 inuh = Math.trunc(inu / 2);
 ir = inu - 2 * inuh;
 arg = (fnu - (inu - ir)) * hpi;
 csgnr = Math.cos(arg);
 csgni = Math.sin(arg);
 if (inuh % 2 === 0) {
 goToLabel = 40;break;
 }
 csgnr = -csgnr;
 csgni = -csgni;
 case 40:

 // -----------------------------------------------------------------------
 // zn is in the right half plane
 // -----------------------------------------------------------------------
 znr = zi;
 zni = -zr;
 if (zi >= 0.0e0) {
 goToLabel = 50;break;
 }
 znr = -znr;
 zni = -zni;
 csgni = -csgni;
 cii = -cii;
 case 50:

 nz = (0, _zbinu.zbinu)(znr, zni, fnu, kode, n, cyr, cyi, rl, fnul, tol, elim, alim);
 if (nz < 0) {
 goToLabel = 130;break;
 }
 nl = n - nz;
 if (nl === 0) break mainExecutionLoop;
 rtol = 1.0e0 / tol;
 ascle = (0, _d1mach.d1mach)(1) * rtol * 1.0e+3;
 for (i = 1; i <= nl; i++) {
 // str = cyr(i)*csgnr - cyi(i)*csgni
 // cyi(i) = cyr(i)*csgni + cyi(i)*csgnr
 // cyr(i) = str
 aa = cyr[i - 1];
 bb = cyi[i - 1];
 atol = 1.0e0;
 if (Math.max(Math.abs(aa), Math.abs(bb)) > ascle) {
 // goToLabel = 55; break;
 } else {
 aa = aa * rtol;
 bb = bb * rtol;
 atol = tol;
 }
 // case 55:

 str = aa * csgnr - bb * csgni;
 sti = aa * csgni + bb * csgnr;
 cyr[i - 1] = str * atol;
 cyi[i - 1] = sti * atol;
 str = -csgni * cii;
 csgni = csgnr * cii;
 csgnr = str;
 }
 break mainExecutionLoop;
 case 130:
 if (nz === -2) {
 goToLabel = 140;break;
 }
 nz = 0;
 ierr = 2;
 break mainExecutionLoop;
 case 140:
 nz = 0;
 ierr = 5;
 break mainExecutionLoop;
 case 260:
 nz = 0;
 ierr = 4;
 default:
 break mainExecutionLoop;
 }
 }
 return [nz, ierr];
} /* eslint-disable no-labels */
/* eslint-disable no-fallthrough */
// **BEGIN PROLOGUE ZBESJ
// **DATE WRITTEN 830501 (YYMMDD)
// **REVISION DATE 890801 (YYMMDD)
// ***PORT TO ECMASCRIPT 201801 (YYYYMM)
// **CATEGORY NO. B5K
// **KEYWORDS J-BESSEL FUNCTION,BESSEL FUNCTION OF COMPLEX ARGUMENT,
// BESSEL FUNCTION OF FIRST KIND
// **AUTHOR (FORTRAN) AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
// **AUTHOR (ECMASCRIPT) ERB, KC, KINGS DISTRIBUTED SYSTEMS
// **PURPOSE TO COMPUTE THE J-BESSEL FUNCTION OF A COMPLEX ARGUMENT
// **DESCRIPTION
//
// ***A DOUBLE PRECISION ROUTINE***
// ON KODE=1, CBESJ COMPUTES AN N MEMBER SEQUENCE OF COMPLEX
// BESSEL FUNCTIONS CY(I)=J(FNU+I-1,Z) FOR REAL, NONNEGATIVE
// ORDERS FNU+I-1, I=1,...,N AND COMPLEX Z IN THE CUT PLANE
// -PI < ARG(Z) <= PI. ON KODE=2, CBESJ RETURNS THE SCALED
// FUNCTIONS
//
// CY(I)=MATH.EXP(-MATH.ABS(Y))*J(FNU+I-1,Z) I = 1,...,N , Y=AIMAG(Z)
//
// WHICH REMOVE THE EXPONENTIAL GROWTH IN BOTH THE UPPER AND
// LOWER HALF PLANES FOR Z TO INFINITY. DEFINITIONS AND NOTATION
// ARE FOUND IN THE NBS HANDBOOK OF MATHEMATICAL FUNCTIONS
// (REF. 1).
//
// INPUT ZR,ZI,FNU ARE DOUBLE PRECISION
// ZR,ZI - Z=CMPLX(ZR,ZI), -PI < ARG(Z) <= PI
// FNU - ORDER OF INITIAL J FUNCTION, FNU >= 0.0E0
// KODE - A PARAMETER TO INDICATE THE SCALING OPTION
// KODE= 1 RETURNS
// CY(I)=J(FNU+I-1,Z), I=1,...,N
// = 2 RETURNS
// CY(I)=J(FNU+I-1,Z)MATH.EXP(-MATH.ABS(Y)), I=1,...,N
// N - NUMBER OF MEMBERS OF THE SEQUENCE, N >= 1
//
// OUTPUT CYR,CYI ARE DOUBLE PRECISION
// CYR,CYI- DOUBLE PRECISION VECTORS WHOSE FIRST N COMPONENTS
// CONTAIN REAL AND IMAGINARY PARTS FOR THE SEQUENCE
// CY(I)=J(FNU+I-1,Z) OR
// CY(I)=J(FNU+I-1,Z)MATH.EXP(-MATH.ABS(Y)) I=1,...,N
// DEPENDING ON KODE, Y=AIMAG(Z).
// NZ - NUMBER OF COMPONENTS SET TO ZERO DUE TO UNDERFLOW,
// NZ= 0 , NORMAL RETURN
// NZ > 0 , LAST NZ COMPONENTS OF CY SET ZERO DUE
// TO UNDERFLOW, CY(I)=CMPLX(0.0E0,0.0E0),
// I = N-NZ+1,...,N
// IERR - ERROR FLAG
// IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
// IERR=1, INPUT ERROR - NO COMPUTATION
// IERR=2, OVERFLOW - NO COMPUTATION, AIMAG(Z)
// TOO LARGE ON KODE=1
// IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE
// BUT LOSSES OF SIGNIFCANCE BY ARGUMENT
// REDUCTION PRODUCE LESS THAN HALF OF MACHINE
// ACCURACY
// IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTA-
// TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-
// CANCE BY ARGUMENT REDUCTION
// IERR=5, ERROR - NO COMPUTATION,
// ALGORITHM TERMINATION CONDITION NOT MET
//
// **LONG DESCRIPTION
//
// THE COMPUTATION IS CARRIED OUT BY THE FORMULA
//
// J(FNU,Z)=MATH.EXP( FNU*PI*MATH.TRUNC(I/2))*I(FNU,-I*Z) AIMAG(Z) >= 0.0
//
// J(FNU,Z)=MATH.EXP(-FNU*PI*MATH.TRUNC(I/2))*I(FNU, I*Z) AIMAG(Z) < 0.0
//
// WHERE I**2 = -1 AND I(FNU,Z) IS THE I BESSEL FUNCTION.
//
// FOR NEGATIVE ORDERS,THE FORMULA
//
// J(-FNU,Z) = J(FNU,Z)*MATH.COS(PI*FNU) - Y(FNU,Z)*MATH.SIN(PI*FNU)
//
// CAN BE USED. HOWEVER,FOR LARGE ORDERS CLOSE TO INTEGERS, THE
// THE FUNCTION CHANGES RADIY(). WHEN FNU IS A LARGE POSITIVE
// INTEGER,THE MAGNITUDE OF J(-FNU,Z)=J(FNU,Z)*MATH.COS(PI*FNU) IS A
// LARGE NEGATIVE POWER OF TEN. BUT WHEN FNU IS NOT AN INTEGER,
// Y(FNU,Z) DOMINATES IN MAGNITUDE WITH A LARGE POSITIVE POWER OF
// TEN AND THE MOST THAT THE SECOND TERM CAN BE REDUCED IS BY
// UNIT ROUNDOFF FROM THE COEFFICIENT. THUS, WIDE CHANGES CAN
// OCCUR WITHIN UNIT ROUNDOFF OF A LARGE INTEGER FOR FNU. HERE,
// LARGE MEANS FNU > CABS(Z).
//
// IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
// MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z OR FNU+N-1 IS
// LARGE, LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR.
// CONSEQUENTLY, IF EITHER ONE EXCEEDS U1=MATH.SQRT(0.5/UR), {
// LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR FLAG
// IERR=3 IS TRIGGERED WHERE UR=MATH.MAX(D1MACH(4),1.0E-18) IS
// DOUBLE PRECISION UNIT ROUNDOFF LIMITED TO 18 DIGITS PRECISION.
// IF EITHER IS LARGER THAN U2=0.5/UR, { ALL SIGNIFICANCE IS
// LOST AND IERR=4. IN ORDER TO USE THE INT FUNCTION, ARGUMENTS
// MUST BE FURTHER RESTRICTED NOT TO EXCEED THE LARGEST MACHINE
// INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF Z AND FNU+N-1 IS
// RESTRICTED BY MATH.MIN(U2,U3). ON 32 BIT MACHINES, U1,U2, AND U3
// ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE PRECISION
// ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE PRECISION
// ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMITING IN
// THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT ONE CAN EXPECT
// TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES, NO DIGITS
// IN SINGLE AND ONLY 7 DIGITS IN DOUBLE PRECISION ARITHMETIC.
// SIMILAR CONSIDERATIONS HOLD FOR OTHER MACHINES.
//
// THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
// BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MATH.MAX(UNIT
// ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
// SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
// ELEMENTARY FUNCTIONS. HERE, S=MATH.MAX(1,MATH.ABS(MATH.LOG10(CABS(Z))),
// MATH.ABS(MATH.LOG10(FNU))) APPROXIMATELY (I.E. S=MATH.MAX(1,MATH.ABS(EXPONENT OF
// CABS(Z),ABS(EXPONENT OF FNU))) ). HOWEVER, THE PHASE ANGLE MAY
// HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
// ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
// SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
// THAN THE OTHER, { ONE CAN EXPECT ONLY MATH.MAX(MATH.ABS(MATH.LOG10(P))-K,
// 0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
// THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
// COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
// BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
// COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
// MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
// THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PMATH.TRUNC(I/2)-P,
// OR -PMATH.TRUNC(I/2)+P.
//
// **REFERENCES HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
// AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
// COMMERCE, 1955.
//
// COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
// BY D. E. AMOS, SAND83-0083, MAY, 1983.
//
// COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
// AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
//
// A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
// ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
// 1018, MAY, 1985
//
// A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
// ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, TRANS.
// MATH. SOFTWARE, 1986
//
// **ROUTINES ED() ZBINU,I1MACH,D1MACH
// **END PROLOGUE ZBESJ
//
// COMPLEX CI,CSGN,CY,Z,ZN
},{"../../utils/fortran-utils/d1mach.js":91,"../../utils/fortran-utils/i1mach.js":92,"./zbinu.js":21}],19:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.zbesk = zbesk;

var _d1mach = require('../../utils/fortran-utils/d1mach.js');

var _i1mach = require('../../utils/fortran-utils/i1mach.js');

var _zabs = require('./zabs.js');

var _zuoik = require('./zuoik.js');

var _zbknu = require('./zbknu.js');

var _zacon = require('./zacon.js');

var _zbunk = require('./zbunk.js');

function zbesk(zr, zi, fnu, kode, n, cyr, cyi) {
 var aa = void 0,
 alim = void 0,
 aln = void 0,
 arg = void 0,
 az = void 0,
 dig = void 0,
 elim = void 0,
 fn = void 0,
 fnul = void 0,
 rl = void 0,
 r1m5 = void 0,
 tol = void 0,
 ufl = void 0,
 bb = void 0,
 ierr = void 0,
 k = void 0,
 k1 = void 0,
 k2 = void 0,
 mr = void 0,
 nn = void 0,
 nuf = void 0,
 nw = void 0,
 nz = void 0;

 var goToLabel = 0;
 mainExecutionLoop: while (true) {
 switch (goToLabel) {
 case 0:
 ierr = 0;
 nz = 0;
 if (zi === 0 && zr === 0) ierr = 1;
 if (fnu < 0.0) ierr = 1;
 if (kode < 1 || kode > 2) ierr = 1;
 if (n < 1) ierr = 1;
 if (ierr !== 0) break mainExecutionLoop;
 nn = n;
 // c-----------------------------------------------------------------------
 // c set parameters related to machine constants.
 // c tol is the approximate unit roundoff limited to 1.0e-18.
 // c elim is the approximate exponential over- and underflow limit.
 // c exp(-elim) < exp(-alim)=exp(-elim)/tol and
 // c exp(elim) > exp(alim)=exp(elim)*tol are intervals near
 // c underflow and overflow limits where scaled arithmetic is done.
 // c rl is the lower boundary of the asymptotic expansion for large z.
 // c dig = number of base 10 digits in tol = 10**(-dig).
 // c fnul is the lower boundary of the asymptotic series for large fnu
 // c-----------------------------------------------------------------------
 tol = Math.max((0, _d1mach.d1mach)(4), 1.0e-18);
 k1 = (0, _i1mach.i1mach)(15);
 k2 = (0, _i1mach.i1mach)(16);
 r1m5 = (0, _d1mach.d1mach)(5);
 k = Math.min(Math.abs(k1), Math.abs(k2));
 elim = 2.303 * k * r1m5 - 3.0;
 k1 = (0, _i1mach.i1mach)(14) - 1;
 aa = r1m5 * k1;
 dig = Math.min(aa, 18.0);
 aa = aa * 2.303;
 alim = elim + Math.max(-aa, -41.45);
 fnul = 10.0 + 6.0 * (dig - 3.0);
 rl = 1.2 * dig + 3.0;
 // c-----------------------------------------------------------------------------
 // c test for proper range
 // c-----------------------------------------------------------------------
 az = (0, _zabs.azabs)(zr, zi);
 fn = fnu + (nn - 1);
 aa = 0.5 / tol;
 bb = (0, _i1mach.i1mach)(9) * 0.5;
 aa = Math.min(aa, bb);
 if (az > aa) {
 goToLabel = 260;break;
 }
 if (fn > aa) {
 goToLabel = 260;break;
 }
 aa = Math.sqrt(aa);
 if (az > aa) ierr = 3;
 if (fn > aa) ierr = 3;
 // c-----------------------------------------------------------------------
 // c overflow test on the last member of the sequence
 // c-----------------------------------------------------------------------
 // c ufl = Math.exp(-elim)
 ufl = (0, _d1mach.d1mach)(1) * 1.0e+3;
 if (az < ufl) {
 goToLabel = 180;break;
 }
 if (fnu > fnul) {
 goToLabel = 80;break;
 }
 if (fn <= 1.0) {
 goToLabel = 60;break;
 }
 if (fn > 2.0) {
 goToLabel = 50;break;
 }
 if (az > tol) {
 goToLabel = 60;break;
 }
 arg = 0.5 * az;
 aln = -fn * Math.log(arg);
 if (aln > elim) {
 goToLabel = 180;break;
 }
 goToLabel = 60;break;
 case 50:
 nuf = (0, _zuoik.zuoik)(zr, zi, fnu, kode, 2, nn, cyr, cyi, tol, elim, alim);
 if (nuf < 0) {
 goToLabel = 180;break;
 }
 nz = nz + nuf;
 nn = nn - nuf;
 // c-----------------------------------------------------------------------
 // c here nn=n or nn=0 since nuf=0,nn, or -1 on return from cuoik
 // c if nuf=nn, then cy(i)=czero for all i
 // c-----------------------------------------------------------------------
 if (nn === 0) {
 goToLabel = 100;break;
 }
 case 60:
 if (zr < 0.0) {
 goToLabel = 70;break;
 }
 // c-----------------------------------------------------------------------
 // c right half plane computation, real(z) >= 0.
 // c-----------------------------------------------------------------------
 nw = (0, _zbknu.zbknu)(zr, zi, fnu, kode, nn, cyr, cyi, tol, elim, alim);
 if (nw < 0) {
 goToLabel = 200;break;
 }
 nz = nw;
 break mainExecutionLoop;
 // c-----------------------------------------------------------------------
 // c left half plane computation
 // c pi/2 < arg(z) <= pi and -pi < arg(z) < -pi/2.
 // c-----------------------------------------------------------------------
 case 70:
 if (nz !== 0) {
 goToLabel = 180;break;
 }
 mr = 1;
 if (zi < 0.0) mr = -1;
 nw = (0, _zacon.zacon)(zr, zi, fnu, kode, mr, nn, cyr, cyi, rl, fnul, tol, elim, alim);
 if (nw < 0) {
 goToLabel = 200;break;
 }
 nz = nw;
 break mainExecutionLoop;
 // c-----------------------------------------------------------------------
 // c uniform asymptotic expansions for fnu > fnul
 // c-----------------------------------------------------------------------
 case 80:
 mr = 0;
 if (zr >= 0.0) {
 goToLabel = 90;break;
 }
 mr = 1;
 if (zi < 0.0) mr = -1;
 case 90:
 nw = (0, _zbunk.zbunk)(zr, zi, fnu, kode, mr, nn, cyr, cyi, tol, elim, alim);
 if (nw < 0) {
 goToLabel = 200;break;
 }
 nz = nz + nw;
 break mainExecutionLoop;
 case 100:
 if (zr < 0.0) {
 goToLabel = 180;break;
 }
 break mainExecutionLoop;
 case 180:
 nz = 0;
 ierr = 2;
 break mainExecutionLoop;
 case 200:
 if (nw === -1) {
 goToLabel = 180;break;
 }
 nz = 0;
 ierr = 5;
 break mainExecutionLoop;
 case 260:
 nz = 0;
 ierr = 4;
 default:
 break mainExecutionLoop;
 }
 }

 return [nz, ierr];
} /* eslint-disable no-labels */
/* eslint-disable no-fallthrough */
// SUBROUTINE ZBESK(ZR, ZI, FNU, KODE, N, CYR, CYI, NZ, IERR)
// ***BEGIN PROLOGUE ZBESK
// ***DATE WRITTEN 830501 (YYMMDD)
// ***REVISION DATE 890801 (YYMMDD)
// ***PORT TO ECMASCRIPT 201801 (YYYYMM)
// ***CATEGORY NO. B5K
// ***KEYWORDS K-BESSEL FUNCTION,COMPLEX BESSEL FUNCTION,
// MODIFIED BESSEL FUNCTION OF THE SECOND KIND,
// BESSEL FUNCTION OF THE THIRD KIND
// ***AUTHOR (FORTRAN) AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
// **AUTHOR (ECMASCRIPT) ERB, KC, KINGS DISTRIBUTED SYSTEMS
// ***PURPOSE TO COMPUTE K-BESSEL FUNCTIONS OF COMPLEX ARGUMENT
// ***DESCRIPTION
//
// ***A DOUBLE PRECISION ROUTINE***
//
// ON KODE=1, CBESK COMPUTES AN N MEMBER SEQUENCE OF COMPLEX
// BESSEL FUNCTIONS CY(J)=K(FNU+J-1,Z) FOR REAL, NONNEGATIVE
// ORDERS FNU+J-1, J=1,...,N AND COMPLEX Z.NE.CMPLX(0.0,0.0)
// IN THE CUT PLANE -PI.LT.ARG(Z).LE.PI. ON KODE=2, CBESK
// RETURNS THE SCALED K FUNCTIONS,
//
// CY(J)=EXP(Z)*K(FNU+J-1,Z) , J=1,...,N,
//
// WHICH REMOVE THE EXPONENTIAL BEHAVIOR IN BOTH THE LEFT AND
// RIGHT HALF PLANES FOR Z TO INFINITY. DEFINITIONS AND
// NOTATION ARE FOUND IN THE NBS HANDBOOK OF MATHEMATICAL
// FUNCTIONS (REF. 1).
//
// INPUT ZR,ZI,FNU ARE DOUBLE PRECISION
// ZR,ZI - Z=CMPLX(ZR,ZI), Z.NE.CMPLX(0.0D0,0.0D0),
// -PI.LT.ARG(Z).LE.PI
// FNU - ORDER OF INITIAL K FUNCTION, FNU.GE.0.0D0
// N - NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
// KODE - A PARAMETER TO INDICATE THE SCALING OPTION
// KODE= 1 RETURNS
// CY(I)=K(FNU+I-1,Z), I=1,...,N
// = 2 RETURNS
// CY(I)=K(FNU+I-1,Z)*EXP(Z), I=1,...,N
//
// OUTPUT CYR,CYI ARE DOUBLE PRECISION
// CYR,CYI- DOUBLE PRECISION VECTORS WHOSE FIRST N COMPONENTS
// CONTAIN REAL AND IMAGINARY PARTS FOR THE SEQUENCE
// CY(I)=K(FNU+I-1,Z), I=1,...,N OR
// CY(I)=K(FNU+I-1,Z)*EXP(Z), I=1,...,N
// DEPENDING ON KODE
// NZ - NUMBER OF COMPONENTS SET TO ZERO DUE TO UNDERFLOW.
// NZ= 0 , NORMAL RETURN
// NZ.GT.0 , FIRST NZ COMPONENTS OF CY SET TO ZERO DUE
// TO UNDERFLOW, CY(I)=CMPLX(0.0D0,0.0D0),
// I=1,...,N WHEN X.GE.0.0. WHEN X.LT.0.0
// NZ STATES ONLY THE NUMBER OF UNDERFLOWS
// IN THE SEQUENCE.
//
// IERR - ERROR FLAG
// IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
// IERR=1, INPUT ERROR - NO COMPUTATION
// IERR=2, OVERFLOW - NO COMPUTATION, FNU IS
// TOO LARGE OR CABS(Z) IS TOO SMALL OR BOTH
// IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE
// BUT LOSSES OF SIGNIFCANCE BY ARGUMENT
// REDUCTION PRODUCE LESS THAN HALF OF MACHINE
// ACCURACY
// IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTA-
// TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-
// CANCE BY ARGUMENT REDUCTION
// IERR=5, ERROR - NO COMPUTATION,
// ALGORITHM TERMINATION CONDITION NOT MET
//
// ***LONG DESCRIPTION
//
// EQUATIONS OF THE REFERENCE ARE IMPLEMENTED FOR SMALL ORDERS
// DNU AND DNU+1.0 IN THE RIGHT HALF PLANE X.GE.0.0. FORWARD
// RECURRENCE GENERATES HIGHER ORDERS. K IS CONTINUED TO THE LEFT
// HALF PLANE BY THE RELATION
//
// K(FNU,Z*EXP(MP)) = EXP(-MP*FNU)*K(FNU,Z)-MP*I(FNU,Z)
// MP=MR*PI*I, MR=+1 OR -1, RE(Z).GT.0, I**2=-1
//
// WHERE I(FNU,Z) IS THE I BESSEL FUNCTION.
//
// FOR LARGE ORDERS, FNU.GT.FNUL, THE K FUNCTION IS COMPUTED
// BY MEANS OF ITS UNIFORM ASYMPTOTIC EXPANSIONS.
//
// FOR NEGATIVE ORDERS, THE FORMULA
//
// K(-FNU,Z) = K(FNU,Z)
//
// CAN BE USED.
//
// CBESK ASSUMES THAT A SIGNIFICANT DIGIT SINH(X) FUNCTION IS
// AVAILABLE.
//
// IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
// MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z OR FNU+N-1 IS
// LARGE, LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR.
// CONSEQUENTLY, IF EITHER ONE EXCEEDS U1=SQRT(0.5/UR), THEN
// LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR FLAG
// IERR=3 IS TRIGGERED WHERE UR=DMAX1(D1MACH(4),1.0D-18) IS
// DOUBLE PRECISION UNIT ROUNDOFF LIMITED TO 18 DIGITS PRECISION.
// IF EITHER IS LARGER THAN U2=0.5/UR, THEN ALL SIGNIFICANCE IS
// LOST AND IERR=4. IN ORDER TO USE THE INT FUNCTION, ARGUMENTS
// MUST BE FURTHER RESTRICTED NOT TO EXCEED THE LARGEST MACHINE
// INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF Z AND FNU+N-1 IS
// RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2, AND U3
// ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE PRECISION
// ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE PRECISION
// ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMITING IN
// THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT ONE CAN EXPECT
// TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES, NO DIGITS
// IN SINGLE AND ONLY 7 DIGITS IN DOUBLE PRECISION ARITHMETIC.
// SIMILAR CONSIDERATIONS HOLD FOR OTHER MACHINES.
//
// THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
// BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
// ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
// SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
// ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
// ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
// CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
// HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
// ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
// SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
// THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
// 0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
// THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
// COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
// BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
// COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
// MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
// THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
// OR -PI/2+P.
//
// ***REFERENCES HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
// AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
// COMMERCE, 1955.
//
// COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
// BY D. E. AMOS, SAND83-0083, MAY, 1983.
//
// COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
// AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983.
//
// A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
// ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
// 1018, MAY, 1985
//
// A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
// ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, TRANS.
// MATH. SOFTWARE, 1986
//
// ***ROUTINES CALLED ZACON,ZBKNU,ZBUNK,ZUOIK,AZABS,I1MACH,D1MACH
// ***END PROLOGUE ZBESK
//
// COMPLEX CY,Z
},{"../../utils/fortran-utils/d1mach.js":91,"../../utils/fortran-utils/i1mach.js":92,"./zabs.js":11,"./zacon.js":13,"./zbknu.js":23,"./zbunk.js":25,"./zuoik.js":44}],20:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});

var _slicedToArray = function () { function sliceIterator(arr, i) { var _arr = []; var _n = true; var _d = false; var _e = undefined; try { for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) { _arr.push(_s.value); if (i && _arr.length === i) break; } } catch (err) { _d = true; _e = err; } finally { try { if (!_n && _i["return"]) _i["return"](); } finally { if (_d) throw _e; } } return _arr; } return function (arr, i) { if (Array.isArray(arr)) { return arr; } else if (Symbol.iterator in Object(arr)) { return sliceIterator(arr, i); } else { throw new TypeError("Invalid attempt to destructure non-iterable instance"); } }; }(); /* eslint-disable no-labels */
/* eslint-disable no-fallthrough */
// **BEGIN PROLOGUE ZBESY
// **DATE WRITTEN 830501 (YYMMDD)
// **REVISION DATE 890801 (YYMMDD)
// ***PORT TO ECMASCRIPT 201801 (YYYYMM)
// **CATEGORY NO. B5K
// **KEYWORDS Y-BESSEL FUNCTION,BESSEL FUNCTION OF COMPLEX ARGUMENT,
// BESSEL FUNCTION OF SECOND KIND
// **AUTHOR (FORTRAN) AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
// **AUTHOR (ECMASCRIPT) ERB, KC, KINGS DISTRIBUTED SYSTEMS
// **PURPOSE TO COMPUTE THE Y-BESSEL FUNCTION OF A COMPLEX ARGUMENT
// **DESCRIPTION
//
// ***A DOUBLE PRECISION ROUTINE***
//
// ON KODE=1, CBESY COMPUTES AN N MEMBER SEQUENCE OF COMPLEX
// BESSEL FUNCTIONS CY(I)=Y(FNU+I-1,Z) FOR REAL, NONNEGATIVE
// ORDERS FNU+I-1, I=1,...,N AND COMPLEX Z IN THE CUT PLANE
// -PI < ARG(Z) <= PI. ON KODE=2, CBESY RETURNS THE SCALED
// FUNCTIONS
//
// CY(I)=MATH.EXP(-MATH.ABS(Y))*Y(FNU+I-1,Z) I = 1,...,N , Y=AIMAG(Z)
//
// WHICH REMOVE THE EXPONENTIAL GROWTH IN BOTH THE UPPER AND
// LOWER HALF PLANES FOR Z TO INFINITY. DEFINITIONS AND NOTATION
// ARE FOUND IN THE NBS HANDBOOK OF MATHEMATICAL FUNCTIONS
// (REF. 1).
//
// INPUT ZR,ZI,FNU ARE DOUBLE PRECISION
// ZR,ZI - Z=CMPLX(ZR,ZI), Z !== CMPLX(0.0E0,0.0E0),
// -PI < ARG(Z) <= PI
// FNU - ORDER OF INITIAL Y FUNCTION, FNU >= 0.0E0
// KODE - A PARAMETER TO INDICATE THE SCALING OPTION
// KODE= 1 RETURNS
// CY(I)=Y(FNU+I-1,Z), I=1,...,N
// = 2 RETURNS
// CY(I)=Y(FNU+I-1,Z)*MATH.EXP(-MATH.ABS(Y)), I=1,...,N
// WHERE Y=AIMAG(Z)
// N - NUMBER OF MEMBERS OF THE SEQUENCE, N >= 1
// CWRKR, - DOUBLE PRECISION WORK VECTORS OF DIMENSION AT
// CWRKI AT LEAST N
//
// OUTPUT CYR,CYI ARE DOUBLE PRECISION
// CYR,CYI- DOUBLE PRECISION VECTORS WHOSE FIRST N COMPONENTS
// CONTAIN REAL AND IMAGINARY PARTS FOR THE SEQUENCE
// CY(I)=Y(FNU+I-1,Z) OR
// CY(I)=Y(FNU+I-1,Z)*MATH.EXP(-MATH.ABS(Y)) I=1,...,N
// DEPENDING ON KODE.
// NZ - NZ=0 , A NORMAL RETURN
// NZ > 0 , NZ COMPONENTS OF CY SET TO ZERO DUE TO
// UNDERFLOW (GENERALLY ON KODE=2)
// IERR - ERROR FLAG
// IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
// IERR=1, INPUT ERROR - NO COMPUTATION
// IERR=2, OVERFLOW - NO COMPUTATION, FNU IS
// TOO LARGE OR CABS(Z) IS TOO SMALL OR BOTH
// IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE
// BUT LOSSES OF SIGNIFCANCE BY ARGUMENT
// REDUCTION PRODUCE LESS THAN HALF OF MACHINE
// ACCURACY
// IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTA-
// TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-
// CANCE BY ARGUMENT REDUCTION
// IERR=5, ERROR - NO COMPUTATION,
// ALGORITHM TERMINATION CONDITION NOT MET
//
// **LONG DESCRIPTION
//
// THE COMPUTATION IS CARRIED OUT BY THE FORMULA
//
// Y(FNU,Z)=0.5*(H(1,FNU,Z)-H(2,FNU,Z))/I
//
// WHERE I**2 = -1 AND THE HANKEL BESSEL FUNCTIONS H(1,FNU,Z)
// AND H(2,FNU,Z) ARE CALCULATED IN CBESH.
//
// FOR NEGATIVE ORDERS,THE FORMULA
//
// Y(-FNU,Z) = Y(FNU,Z)*MATH.COS(PI*FNU) + J(FNU,Z)*MATH.SIN(PI*FNU)
//
// CAN BE USED. HOWEVER,FOR LARGE ORDERS CLOSE TO HALF ODD
// INTEGERS THE FUNCTION CHANGES RADIY(). WHEN FNU IS A LARGE
// POSITIVE HALF ODD INTEGER,THE MAGNITUDE OF Y(-FNU,Z)=J(FNU,Z)*
// MATH.SIN(PI*FNU) IS A LARGE NEGATIVE POWER OF TEN. BUT WHEN FNU IS
// NOT A HALF ODD INTEGER, Y(FNU,Z) DOMINATES IN MAGNITUDE WITH A
// LARGE POSITIVE POWER OF TEN AND THE MOST THAT THE SECOND TERM
// CAN BE REDUCED IS BY UNIT ROUNDOFF FROM THE COEFFICIENT. THUS,
// WIDE CHANGES CAN OCCUR WITHIN UNIT ROUNDOFF OF A LARGE HALF
// ODD INTEGER. HERE, LARGE MEANS FNU > CABS(Z).
//
// IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
// MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z OR FNU+N-1 IS
// LARGE, LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR.
// CONSEQUENTLY, IF EITHER ONE EXCEEDS U1=MATH.SQRT(0.5/UR), {
// LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR FLAG
// IERR=3 IS TRIGGERED WHERE UR=MATH.MAX(D1MACH(4),1.0E-18) IS
// DOUBLE PRECISION UNIT ROUNDOFF LIMITED TO 18 DIGITS PRECISION.
// IF EITHER IS LARGER THAN U2=0.5/UR, { ALL SIGNIFICANCE IS
// LOST AND IERR=4. IN ORDER TO USE THE INT FUNCTION, ARGUMENTS
// MUST BE FURTHER RESTRICTED NOT TO EXCEED THE LARGEST MACHINE
// INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF Z AND FNU+N-1 IS
// RESTRICTED BY MATH.MIN(U2,U3). ON 32 BIT MACHINES, U1,U2, AND U3
// ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE PRECISION
// ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE PRECISION
// ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMITING IN
// THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT ONE CAN EXPECT
// TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES, NO DIGITS
// IN SINGLE AND ONLY 7 DIGITS IN DOUBLE PRECISION ARITHMETIC.
// SIMILAR CONSIDERATIONS HOLD FOR OTHER MACHINES.
//
// THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
// BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MATH.MAX(UNIT
// ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
// SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
// ELEMENTARY FUNCTIONS. HERE, S=MATH.MAX(1,MATH.ABS(MATH.LOG10(CABS(Z))),
// MATH.ABS(MATH.LOG10(FNU))) APPROXIMATELY (I.E. S=MATH.MAX(1,MATH.ABS(EXPONENT OF
// CABS(Z),ABS(EXPONENT OF FNU)) )). HOWEVER, THE PHASE ANGLE MAY
// HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
// ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
// SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
// THAN THE OTHER, { ONE CAN EXPECT ONLY MATH.MAX(MATH.ABS(MATH.LOG10(P))-K,
// 0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
// THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
// COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
// BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
// COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
// MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
// THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PMATH.TRUNC(I/2)-P,
// OR -PMATH.TRUNC(I/2)+P.
//
// **REFERENCES HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
// AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
// COMMERCE, 1955.
//
// COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
// BY D. E. AMOS, SAND83-0083, MAY, 1983.
//
// COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
// AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
//
// A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
// ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
// 1018, MAY, 1985
//
// A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
// ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, TRANS.
// MATH. SOFTWARE, 1986
//
// **ROUTINES CALLED ZBESH,I1MACH,D1MACH
// **END PROLOGUE ZBESY
//
// COMPLEX CWRK,CY,C1,C2,EX,HCI,Z,ZU,ZV


exports.zbesy = zbesy;

var _d1mach = require('../../utils/fortran-utils/d1mach.js');

var _i1mach = require('../../utils/fortran-utils/i1mach.js');

var _zbesh5 = require('./zbesh.js');

function zbesy(zr, zi, fnu, kode, n, cyr, cyi) {
 var c1i = void 0,
 c1r = void 0,
 c2i = void 0,
 c2r = void 0,
 elim = void 0,
 exi = void 0,
 exr = void 0,
 ey = void 0,
 hcii = void 0,
 sti = void 0,
 str = void 0,
 tay = void 0,
 ascle = void 0,
 rtol = void 0,
 atol = void 0,
 aa = void 0,
 bb = void 0,
 tol = void 0,
 i = void 0,
 k = void 0,
 k1 = void 0,
 k2 = void 0,
 nz1 = void 0,
 nz2 = void 0,
 cwrkr = void 0,
 cwrki = void 0,
 ierr = void 0,
 nz = void 0,
 r1m5 = void 0;

 cwrkr = new Array(n);
 cwrki = new Array(n);

 var goToLabel = 0;
 mainExecutionLoop: while (true) {
 switch (goToLabel) {
 case 0:
 ierr = 0;
 nz = 0;
 if (zr === 0.0e0 && zi === 0.0e0) ierr = 1;
 if (fnu < 0.0e0) ierr = 1;
 if (kode < 1 || kode > 2) ierr = 1;
 if (n < 1) ierr = 1;
 if (ierr !== 0) break mainExecutionLoop;
 hcii = 0.5e0;

 var _zbesh = (0, _zbesh5.zbesh)(zr, zi, fnu, kode, 1, n, cyr, cyi);

 var _zbesh2 = _slicedToArray(_zbesh, 2);

 nz1 = _zbesh2[0];
 ierr = _zbesh2[1];

 if (ierr !== 0 && ierr !== 3) {
 goToLabel = 170;break;
 }

 var _zbesh3 = (0, _zbesh5.zbesh)(zr, zi, fnu, kode, 2, n, cwrkr, cwrki);

 var _zbesh4 = _slicedToArray(_zbesh3, 2);

 nz2 = _zbesh4[0];
 ierr = _zbesh4[1];

 if (ierr !== 0 && ierr !== 3) {
 goToLabel = 170;break;
 }
 nz = Math.min(nz1, nz2);
 if (kode === 2) {
 goToLabel = 60;break;
 }
 for (i = 1; i <= n; i++) {
 str = cwrkr[i - 1] - cyr[i - 1];
 sti = cwrki[i - 1] - cyi[i - 1];
 cyr[i - 1] = -sti * hcii;
 cyi[i - 1] = str * hcii;
 }
 break mainExecutionLoop;
 case 60:

 tol = Math.max((0, _d1mach.d1mach)(4), 1.0e-18);
 k1 = (0, _i1mach.i1mach)(15);
 k2 = (0, _i1mach.i1mach)(16);
 k = Math.min(Math.abs(k1), Math.abs(k2));
 r1m5 = (0, _d1mach.d1mach)(5);
 // -----------------------------------------------------------------------
 // elim is the approximate exponential under- and overflow limit
 // -----------------------------------------------------------------------
 elim = 2.303e0 * (k * r1m5 - 3.0e0);
 exr = Math.cos(zr);
 exi = Math.sin(zr);
 ey = 0.0e0;
 tay = Math.abs(zi + zi);
 if (tay < elim) ey = Math.exp(-tay);
 if (zi < 0.0e0) {
 goToLabel = 90;break;
 }
 c1r = exr * ey;
 c1i = exi * ey;
 c2r = exr;
 c2i = -exi;
 case 70:

 nz = 0;
 rtol = 1.0e0 / tol;
 ascle = (0, _d1mach.d1mach)(1) * rtol * 1.0e+3;
 for (i = 1; i <= n; i++) {
 // str = c1r*cyr(i) - c1i*cyi(i)
 // sti = c1r*cyi(i) + c1i*cyr(i)
 // str = -str + c2r*cwrkr(i) - c2i*cwrki(i)
 // sti = -sti + c2r*cwrki(i) + c2i*cwrkr(i)
 // cyr(i) = -sti*hcii
 // cyi(i) = str*hcii
 aa = cwrkr[i - 1];
 bb = cwrki[i - 1];
 atol = 1.0e0;
 if (Math.max(Math.abs(aa), Math.abs(bb)) > ascle) {
 // goToLabel = 75; break;
 } else {
 aa = aa * rtol;
 bb = bb * rtol;
 atol = tol;
 }
 // case 75:

 str = (aa * c2r - bb * c2i) * atol;
 sti = (aa * c2i + bb * c2r) * atol;
 aa = cyr[i - 1];
 bb = cyi[i - 1];
 atol = 1.0e0;
 if (Math.max(Math.abs(aa), Math.abs(bb)) > ascle) {
 // goToLabel = 85; break;
 } else {
 aa = aa * rtol;
 bb = bb * rtol;
 atol = tol;
 }
 // case 85:

 str = str - (aa * c1r - bb * c1i) * atol;
 sti = sti - (aa * c1i + bb * c1r) * atol;
 cyr[i - 1] = -sti * hcii;
 cyi[i - 1] = str * hcii;
 if (str === 0.0e0 && sti === 0.0e0 && ey === 0.0e0) {
 nz = nz + 1;
 }
 }
 break mainExecutionLoop;
 case 90:

 c1r = exr;
 c1i = exi;
 c2r = exr * ey;
 c2i = -exi * ey;
 goToLabel = 70;break;
 case 170:

 nz = 0;
 default:
 break mainExecutionLoop;
 }
 }
 return [nz, ierr];
}
},{"../../utils/fortran-utils/d1mach.js":91,"../../utils/fortran-utils/i1mach.js":92,"./zbesh.js":16}],21:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});

var _slicedToArray = function () { function sliceIterator(arr, i) { var _arr = []; var _n = true; var _d = false; var _e = undefined; try { for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) { _arr.push(_s.value); if (i && _arr.length === i) break; } } catch (err) { _d = true; _e = err; } finally { try { if (!_n && _i["return"]) _i["return"](); } finally { if (_d) throw _e; } } return _arr; } return function (arr, i) { if (Array.isArray(arr)) { return arr; } else if (Symbol.iterator in Object(arr)) { return sliceIterator(arr, i); } else { throw new TypeError("Invalid attempt to destructure non-iterable instance"); } }; }(); /* eslint-disable no-labels */
/* eslint-disable no-fallthrough */
// SUBROUTINE ZBINU(ZR, ZI, FNU, KODE, N, CYR, CYI, NZ, RL, FNUL,
// * TOL, ELIM, ALIM)
// ***BEGIN PROLOGUE ZBINU
// ***REFER TO ZBESH,ZBESI,ZBESJ,ZBESK,ZAIRY,ZBIRY
//
// ZBINU COMPUTES THE I FUNCTION IN THE RIGHT HALF Z PLANE
//
// ***ROUTINES CALLED AZABS,ZASYI,ZBUNI,ZMLRI,ZSERI,ZUOIK,ZWRSK
// ***END PROLOGUE ZBINU


exports.zbinu = zbinu;

var _zabs = require('./zabs.js');

var _zasyi = require('./zasyi.js');

var _zbuni3 = require('./zbuni.js');

var _zmlri = require('./zmlri.js');

var _zseri = require('./zseri.js');

var _zuoik = require('./zuoik.js');

var _zwrsk = require('./zwrsk.js');

function zbinu(zr, zi, fnu, kode, n, cyr, cyi, rl, fnul, tol, elim, alim) {
 var az = void 0,
 cwi = void 0,
 cwr = void 0,
 dfnu = void 0,
 zeroi = void 0,
 zeror = void 0,
 i = void 0,
 inw = void 0,
 nlast = void 0,
 nn = void 0,
 nui = void 0,
 nw = void 0,
 nz = void 0;
 cwr = new Array(2);
 cwi = new Array(2);
 zeror = 0;
 zeroi = 0;


 var goToLabel = 0;
 mainExecutionLoop: while (true) {
 switch (goToLabel) {
 case 0:
 nz = 0;
 az = (0, _zabs.azabs)(zr, zi);
 nn = n;
 dfnu = fnu + (n - 1);
 if (az <= 2.0) {
 goToLabel = 10;break;
 }
 if (az * az * 0.25 > dfnu + 1.0) {
 goToLabel = 20;break;
 }
 case 10:
 // c-----------------------------------------------------------------------
 // c power series
 // c-----------------------------------------------------------------------
 nw = (0, _zseri.zseri)(zr, zi, fnu, kode, nn, cyr, cyi, tol, elim, alim);
 inw = Math.abs(nw);
 nz = nz + inw;
 nn = nn - inw;
 if (nn === 0) break mainExecutionLoop;
 if (nw >= 0) {
 goToLabel = 120;break;
 }
 dfnu = fnu + (nn - 1);
 case 20:
 if (az < rl) {
 goToLabel = 40;break;
 }
 if (dfnu <= 1.0) {
 goToLabel = 30;break;
 }
 if (az + az < dfnu * dfnu) {
 goToLabel = 50;break;
 }
 // c-----------------------------------------------------------------------
 // c asymptotic expansion for large z
 // c-----------------------------------------------------------------------
 case 30:
 nw = (0, _zasyi.zasyi)(zr, zi, fnu, kode, nn, cyr, cyi, rl, tol, elim, alim);
 if (nw < 0) {
 goToLabel = 130;break;
 }
 goToLabel = 120;break;
 case 40:
 if (dfnu <= 1.0) {
 goToLabel = 70;break;
 }
 case 50:
 // c-----------------------------------------------------------------------
 // c overflow and underflow test on i sequence for miller algorithm
 // c-----------------------------------------------------------------------
 nw = (0, _zuoik.zuoik)(zr, zi, fnu, kode, 1, nn, cyr, cyi, tol, elim, alim);
 if (nw < 0) {
 goToLabel = 130;break;
 }
 nz = nz + nw;
 nn = nn - nw;
 if (nn === 0) break mainExecutionLoop;
 dfnu = fnu + (nn - 1);
 if (dfnu > fnul) {
 goToLabel = 110;break;
 }
 if (az > fnul) {
 goToLabel = 110;break;
 }
 case 60:
 if (az > rl) {
 goToLabel = 80;break;
 }
 case 70:
 // c-----------------------------------------------------------------------
 // c miller algorithm normalized by the series
 // c-----------------------------------------------------------------------
 nw = (0, _zmlri.zmlri)(zr, zi, fnu, kode, nn, cyr, cyi, tol);
 if (nw < 0) {
 goToLabel = 130;break;
 }
 goToLabel = 120;break;
 case 80:
 // c-----------------------------------------------------------------------
 // c miller algorithm normalized by the wronskian
 // c-----------------------------------------------------------------------
 // c-----------------------------------------------------------------------
 // c overflow test on k functions used in wronskian
 // c-----------------------------------------------------------------------
 nw = (0, _zuoik.zuoik)(zr, zi, fnu, kode, 2, 2, cwr, cwi, tol, elim, alim);
 if (nw >= 0) {
 goToLabel = 100;break;
 }
 nz = nn;
 // do 90 i=1,nn
 for (i = 1; i <= nn; i++) {
 cyr[i - 1] = zeror;
 cyi[i - 1] = zeroi;
 }
 // 90 continue
 break mainExecutionLoop;
 case 100:
 if (nw > 0) {
 goToLabel = 130;break;
 }
 nw = (0, _zwrsk.zwrsk)(zr, zi, fnu, kode, nn, cyr, cyi, cwr, cwi, tol, elim, alim);
 if (nw < 0) {
 goToLabel = 130;break;
 }
 goToLabel = 120;break;
 case 110:
 // c-----------------------------------------------------------------------
 // c increment fnu+nn-1 up to fnul, compute and recur backward
 // c-----------------------------------------------------------------------
 nui = Math.trunc(fnul - dfnu) + 1;
 nui = Math.max(nui, 0);

 var _zbuni = (0, _zbuni3.zbuni)(zr, zi, fnu, kode, nn, cyr, cyi, nui, nlast, fnul, tol, elim, alim);

 var _zbuni2 = _slicedToArray(_zbuni, 2);

 nw = _zbuni2[0];
 nlast = _zbuni2[1];

 if (nw < 0) {
 goToLabel = 130;break;
 }
 nz = nz + nw;
 if (nlast === 0) {
 goToLabel = 120;break;
 }
 nn = nlast;
 goToLabel = 60;break;
 case 120:
 break mainExecutionLoop;
 case 130:
 nz = -1;
 if (nw === -2) nz = -2;
 default:
 break mainExecutionLoop;
 }
 }
 return nz;
}
},{"./zabs.js":11,"./zasyi.js":15,"./zbuni.js":24,"./zmlri.js":30,"./zseri.js":34,"./zuoik.js":44,"./zwrsk.js":45}],22:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});

var _slicedToArray = function () { function sliceIterator(arr, i) { var _arr = []; var _n = true; var _d = false; var _e = undefined; try { for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) { _arr.push(_s.value); if (i && _arr.length === i) break; } } catch (err) { _d = true; _e = err; } finally { try { if (!_n && _i["return"]) _i["return"](); } finally { if (_d) throw _e; } } return _arr; } return function (arr, i) { if (Array.isArray(arr)) { return arr; } else if (Symbol.iterator in Object(arr)) { return sliceIterator(arr, i); } else { throw new TypeError("Invalid attempt to destructure non-iterable instance"); } }; }(); /* eslint-disable no-labels */
/* eslint-disable no-fallthrough */
// **BEGIN PROLOGUE ZBIRY
// **DATE WRITTEN 830501 (YYMMDD)
// **REVISION DATE 890801 (YYMMDD)
// ***PORT TO ECMASCRIPT 201801 (YYYYMM)
// **CATEGORY NO. B5K
// **KEYWORDS AIRY FUNCTION, BESSEL FUNCTIONS OF ORDER ONE THIRD
// **AUTHOR (FORTRAN) AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
// **AUTHOR (ECMASCRIPT) ERB, KC, KINGS DISTRIBUTED SYSTEMS
// **PURPOSE TO COMPUTE AIRY FUNCTIONS BI(Z) AND DBI(Z) FOR COMPLEX Z
// **DESCRIPTION
//
// ***A DOUBLE PRECISION ROUTINE***
// ON KODE=1, CBIRY COMPUTES THE COMPLEX AIRY FUNCTION BI(Z) OR
// ITS DERIVATIVE DBI(Z)/DZ ON ID=0 OR ID=1 RESPECTIVELY. ON
// KODE=2, A SCALING OPTION CEXP(-AXZTA)*BI(Z) OR CEXP(-AXZTA)*
// DBI(Z)/DZ IS PROVIDED TO REMOVE THE EXPONENTIAL BEHAVIOR IN
// BOTH THE LEFT AND RIGHT HALF PLANES WHERE
// ZTA=(MATH.TRUNC(2/3))*Z*CSQRT(Z)=CMPLX(XZTA,YZTA) AND AXZTA=MATH.ABS(XZTA).
// DEFINITIONS AND NOTATION ARE FOUND IN THE NBS HANDBOOK OF
// MATHEMATICAL FUNCTIONS (REF. 1).
//
// INPUT ZR,ZI ARE DOUBLE PRECISION
// ZR,ZI - Z=CMPLX(ZR,ZI)
// ID - ORDER OF DERIVATIVE, ID=0 OR ID=1
// KODE - A PARAMETER TO INDICATE THE SCALING OPTION
// KODE= 1 RETURNS
// BI=BI(Z) ON ID=0 OR
// BI=DBI(Z)/DZ ON ID=1
// = 2 RETURNS
// BI=CEXP(-AXZTA)*BI(Z) ON ID=0 OR
// BI=CEXP(-AXZTA)*DBI(Z)/DZ ON ID=1 WHERE
// ZTA=(MATH.TRUNC(2/3))*Z*CSQRT(Z)=CMPLX(XZTA,YZTA)
// AND AXZTA=MATH.ABS(XZTA)
//
// OUTPUT BIR,BII ARE DOUBLE PRECISION
// BIR,BII- COMPLEX ANSWER DEPENDING ON THE CHOICES FOR ID AND
// KODE
// IERR - ERROR FLAG
// IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
// IERR=1, INPUT ERROR - NO COMPUTATION
// IERR=2, OVERFLOW - NO COMPUTATION, REAL(Z)
// TOO LARGE ON KODE=1
// IERR=3, CABS(Z) LARGE - COMPUTATION COMPLETED
// LOSSES OF SIGNIFCANCE BY ARGUMENT REDUCTION
// PRODUCE LESS THAN HALF OF MACHINE ACCURACY
// IERR=4, CABS(Z) TOO LARGE - NO COMPUTATION
// COMPLETE LOSS OF ACCURACY BY ARGUMENT
// REDUCTION
// IERR=5, ERROR - NO COMPUTATION,
// ALGORITHM TERMINATION CONDITION NOT MET
//
// **LONG DESCRIPTION
//
// BI AND DBI ARE COMPUTED FOR CABS(Z) > 1.0 FROM THE I BESSEL
// FUNCTIONS BY
//
// BI(Z)=C*MATH.SQRT(Z)*( I(-MATH.TRUNC(1/3),ZTA) + I(MATH.TRUNC(1/3),ZTA) )
// DBI(Z)=C * Z * ( I(-MATH.TRUNC(2/3),ZTA) + I(MATH.TRUNC(2/3),ZTA) )
// C=1.0/MATH.SQRT(3.0)
// ZTA=(MATH.TRUNC(2/3))*Z**(MATH.TRUNC(3/2))
//
// WITH THE POWER SERIES FOR CABS(Z) <= 1.0.
//
// IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
// MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z IS LARGE, LOSSES
// OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR. CONSEQUENTLY, IF
// THE MAGNITUDE OF ZETA=(MATH.TRUNC(2/3))*Z**1.5 EXCEEDS U1=MATH.SQRT(0.5/UR),
// { LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR
// FLAG IERR=3 IS TRIGGERED WHERE UR=MATH.MAX(D1MACH(4),1.0E-18) IS
// DOUBLE PRECISION UNIT ROUNDOFF LIMITED TO 18 DIGITS PRECISION.
// ALSO, IF THE MAGNITUDE OF ZETA IS LARGER THAN U2=0.5/UR, {
// ALL SIGNIFICANCE IS LOST AND IERR=4. IN ORDER TO USE THE INT
// FUNCTION, ZETA MUST BE FURTHER RESTRICTED NOT TO EXCEED THE
// LARGEST INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF ZETA
// MUST BE RESTRICTED BY MATH.MIN(U2,U3). ON 32 BIT MACHINES, U1,U2,
// AND U3 ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE
// PRECISION ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE
// PRECISION ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMIT-
// ING IN THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT THE MAG-
// NITUDE OF Z CANNOT EXCEED 3.1E+4 IN SINGLE AND 2.1E+6 IN
// DOUBLE PRECISION ARITHMETIC. THIS ALSO MEANS THAT ONE CAN
// EXPECT TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES,
// NO DIGITS IN SINGLE PRECISION AND ONLY 7 DIGITS IN DOUBLE
// PRECISION ARITHMETIC. SIMILAR CONSIDERATIONS HOLD FOR OTHER
// MACHINES.
//
// THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
// BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MATH.MAX(UNIT
// ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
// SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
// ELEMENTARY FUNCTIONS. HERE, S=MATH.MAX(1,MATH.ABS(MATH.LOG10(CABS(Z))),
// MATH.ABS(MATH.LOG10(FNU))) APPROXIMATELY (I.E. S=MATH.MAX(1,MATH.ABS(EXPONENT OF
// CABS(Z),ABS(EXPONENT OF FNU))) ). HOWEVER, THE PHASE ANGLE MAY
// HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
// ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
// SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
// THAN THE OTHER, { ONE CAN EXPECT ONLY MATH.MAX(MATH.ABS(MATH.LOG10(P))-K,
// 0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
// THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
// COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
// BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
// COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
// MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
// THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
// OR -PI/2+P.
//
// **REFERENCES HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
// AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
// COMMERCE, 1955.
//
// COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
// AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
//
// A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
// ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
// 1018, MAY, 1985
//
// A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
// ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, TRANS.
// MATH. SOFTWARE, 1986
//
// **ROUTINES ED() ZBINU,AZABS,ZDIV,AZSQRT,D1MACH,I1MACH
// **END PROLOGUE ZBIRY
// COMPLEX BI,CONE,CSQ,CY,S1,S2,TRM1,TRM2,Z,ZTA,Z3


exports.zbiry = zbiry;

var _i1mach = require('../../utils/fortran-utils/i1mach.js');

var _d1mach = require('../../utils/fortran-utils/d1mach.js');

var _zbinu = require('./zbinu.js');

var _zdiv3 = require('./zdiv.js');

var _zsqrt = require('./zsqrt.js');

var _zabs = require('./zabs.js');

function zbiry(zr, zi, id, kode) {
 var aa = void 0,
 ad = void 0,
 ak = void 0,
 alim = void 0,
 atrm = void 0,
 az = void 0,
 az3 = void 0,
 bb = void 0,
 bk = void 0,
 bii = void 0,
 bir = void 0,
 cc = void 0,
 ck = void 0,
 coef = void 0,
 conei = void 0,
 coner = void 0,
 csqi = void 0,
 csqr = void 0,
 cyi = void 0,
 cyr = void 0,
 c1 = void 0,
 c2 = void 0,
 dig = void 0,
 dk = void 0,
 d1 = void 0,
 d2 = void 0,
 eaa = void 0,
 elim = void 0,
 fid = void 0,
 fmr = void 0,
 fnu = void 0,
 fnul = void 0,
 pi = void 0,
 rl = void 0,
 r1m5 = void 0,
 sfac = void 0,
 sti = void 0,
 str = void 0,
 s1i = void 0,
 s1r = void 0,
 s2i = void 0,
 s2r = void 0,
 tol = void 0,
 trm1i = void 0,
 trm1r = void 0,
 trm2i = void 0,
 trm2r = void 0,
 tth = void 0,
 ztai = void 0,
 ztar = void 0,
 z3i = void 0,
 z3r = void 0,
 k = void 0,
 k1 = void 0,
 k2 = void 0,
 nz = void 0,
 ierr = void 0;

 cyr = new Float64Array(2);
 cyi = new Float64Array(2);
 var goToLabel = 0;
 mainExecutionLoop: while (true) {
 switch (goToLabel) {
 case 0:
 tth = 6.66666666666666667e-01;
 c1 = 6.14926627446000736e-01;
 c2 = 4.48288357353826359e-01;
 coef = 5.77350269189625765e-01;
 pi = 3.14159265358979324e+00;

 //* **first executable statement zbiry
 coner = 1.0e0;
 conei = 0.0e0;
 ierr = 0;
 nz = 0;
 if (id < 0 || id > 1) ierr = 1;
 if (kode < 1 || kode > 2) ierr = 1;
 if (ierr !== 0) return;
 az = (0, _zabs.azabs)(zr, zi);
 tol = Math.max((0, _d1mach.d1mach)(4), 1.0e-18);
 fid = id;
 if (az > 1.0e0) {
 goToLabel = 70;break;
 }
 // -----------------------------------------------------------------------
 // power series for cabs(z) <= 1.
 // -----------------------------------------------------------------------
 s1r = coner;
 s1i = conei;
 s2r = coner;
 s2i = conei;
 if (az < tol) {
 goToLabel = 130;break;
 }
 aa = az * az;
 if (aa < tol / az) {
 goToLabel = 40;break;
 }
 trm1r = coner;
 trm1i = conei;
 trm2r = coner;
 trm2i = conei;
 atrm = 1.0e0;
 str = zr * zr - zi * zi;
 sti = zr * zi + zi * zr;
 z3r = str * zr - sti * zi;
 z3i = str * zi + sti * zr;
 az3 = az * aa;
 ak = 2.0e0 + fid;
 bk = 3.0e0 - fid - fid;
 ck = 4.0e0 - fid;
 dk = 3.0e0 + fid + fid;
 d1 = ak * dk;
 d2 = bk * ck;
 ad = Math.min(d1, d2);
 ak = 24.0e0 + 9.0e0 * fid;
 bk = 30.0e0 - 9.0e0 * fid;
 for (k = 1; k <= 25; k++) {
 str = (trm1r * z3r - trm1i * z3i) / d1;
 trm1i = (trm1r * z3i + trm1i * z3r) / d1;
 trm1r = str;
 s1r = s1r + trm1r;
 s1i = s1i + trm1i;
 str = (trm2r * z3r - trm2i * z3i) / d2;
 trm2i = (trm2r * z3i + trm2i * z3r) / d2;
 trm2r = str;
 s2r = s2r + trm2r;
 s2i = s2i + trm2i;
 atrm = atrm * az3 / ad;
 d1 = d1 + ak;
 d2 = d2 + bk;
 ad = Math.min(d1, d2);
 if (atrm < tol * ad) break;
 ak = ak + 18.0e0;
 bk = bk + 18.0e0;
 }
 case 40:

 if (id === 1) {
 goToLabel = 50;break;
 }
 bir = c1 * s1r + c2 * (zr * s2r - zi * s2i);
 bii = c1 * s1i + c2 * (zr * s2i + zi * s2r);
 if (kode === 1) break mainExecutionLoop;

 var _azsqrt = (0, _zsqrt.azsqrt)(zr, zi);

 var _azsqrt2 = _slicedToArray(_azsqrt, 2);

 str = _azsqrt2[0];
 sti = _azsqrt2[1];

 ztar = tth * (zr * str - zi * sti);
 ztai = tth * (zr * sti + zi * str);
 aa = ztar;
 aa = -Math.abs(aa);
 eaa = Math.exp(aa);
 bir = bir * eaa;
 bii = bii * eaa;
 return;
 case 50:

 bir = s2r * c2;
 bii = s2i * c2;
 if (az <= tol) {
 goToLabel = 60;break;
 }
 cc = c1 / (1.0e0 + fid);
 str = s1r * zr - s1i * zi;
 sti = s1r * zi + s1i * zr;
 bir = bir + cc * (str * zr - sti * zi);
 bii = bii + cc * (str * zi + sti * zr);
 case 60:

 if (kode === 1) return;

 var _azsqrt3 = (0, _zsqrt.azsqrt)(zr, zi);

 var _azsqrt4 = _slicedToArray(_azsqrt3, 2);

 str = _azsqrt4[0];
 sti = _azsqrt4[1];

 ztar = tth * (zr * str - zi * sti);
 ztai = tth * (zr * sti + zi * str);
 aa = ztar;
 aa = -Math.abs(aa);
 eaa = Math.exp(aa);
 bir = bir * eaa;
 bii = bii * eaa;
 return;
 // -----------------------------------------------------------------------
 // case for cabs(z) > 1.0
 // -----------------------------------------------------------------------
 case 70:

 fnu = (1.0e0 + fid) / 3.0e0;
 // -----------------------------------------------------------------------
 // set parameters related to machine constants.
 // tol is the approximate unit roundoff limited to 1.0e-18.
 // elim is the approximate exponential over- and underflow limit.
 // Math.exp(-elim) < Math.exp(-alim)=Math.exp(-elim)/tol and
 // Math.exp(elim) > Math.exp(alim)=Math.exp(elim)*tol are intervals near
 // underflow and overflow limits where scaled arithmetic is done.
 // rl is the lower boundary of the asymptotic expansion for large z.
 // dig = number of base 10 digits in tol = 10**(-dig).
 // fnul is the lower boundary of the asymptotic series for large fnu.
 // -----------------------------------------------------------------------
 k1 = (0, _i1mach.i1mach)(15);
 k2 = (0, _i1mach.i1mach)(16);
 r1m5 = (0, _d1mach.d1mach)(5);
 k = Math.min(Math.abs(k1), Math.abs(k2));
 elim = 2.303e0 * (k * r1m5 - 3.0e0);
 k1 = (0, _i1mach.i1mach)(14) - 1;
 aa = r1m5 * k1;
 dig = Math.min(aa, 18.0e0);
 aa = aa * 2.303e0;
 alim = elim + Math.max(-aa, -41.45e0);
 rl = 1.2e0 * dig + 3.0e0;
 fnul = 10.0e0 + 6.0e0 * (dig - 3.0e0);
 // -----------------------------------------------------------------------
 // test for range
 // -----------------------------------------------------------------------
 aa = 0.5e0 / tol;
 bb = (0, _i1mach.i1mach)(9) * 0.5e0;
 aa = Math.min(aa, bb);
 aa = aa ** tth;
 if (az > aa) {
 goToLabel = 260;break;
 }
 aa = Math.sqrt(aa);
 if (az > aa) ierr = 3;

 var _azsqrt5 = (0, _zsqrt.azsqrt)(zr, zi);

 var _azsqrt6 = _slicedToArray(_azsqrt5, 2);

 csqr = _azsqrt6[0];
 csqi = _azsqrt6[1];

 ztar = tth * (zr * csqr - zi * csqi);
 ztai = tth * (zr * csqi + zi * csqr);
 // -----------------------------------------------------------------------
 // re(zta) <= 0 when re(z) < 0, especially when im(z) is small
 // -----------------------------------------------------------------------
 sfac = 1.0e0;
 ak = ztai;
 if (zr >= 0.0e0) {
 goToLabel = 80;break;
 }
 bk = ztar;
 ck = -Math.abs(bk);
 ztar = ck;
 ztai = ak;
 case 80:

 if (zi !== 0.0e0 || zr > 0.0e0) {
 goToLabel = 90;break;
 }
 ztar = 0.0e0;
 ztai = ak;
 case 90:

 aa = ztar;
 if (kode === 2) {
 goToLabel = 100;break;
 }
 // -----------------------------------------------------------------------
 // overflow test
 // -----------------------------------------------------------------------
 bb = Math.abs(aa);
 if (bb < alim) {
 goToLabel = 100;break;
 }
 bb = bb + 0.25e0 * Math.log(az);
 sfac = tol;
 if (bb > elim) {
 goToLabel = 190;break;
 }
 case 100:

 fmr = 0.0e0;
 if (aa >= 0.0e0 && zr > 0.0e0) {
 goToLabel = 110;break;
 }
 fmr = pi;
 if (zi < 0.0e0) fmr = -pi;
 ztar = -ztar;
 ztai = -ztai;
 case 110:

 // -----------------------------------------------------------------------
 // aa=factor for analytic continuation of i(fnu,zta)
 // kode=2 returns Math.exp(-Math.abs(xzta))*i(fnu,zta) from cbesi
 // -----------------------------------------------------------------------
 nz = (0, _zbinu.zbinu)(ztar, ztai, fnu, kode, 1, cyr, cyi, rl, fnul, tol, elim, alim);
 if (nz < 0) {
 goToLabel = 200;break;
 }
 aa = fmr * fnu;
 z3r = sfac;
 str = Math.cos(aa);
 sti = Math.sin(aa);
 s1r = (str * cyr[0] - sti * cyi[0]) * z3r;
 s1i = (str * cyi[0] + sti * cyr[0]) * z3r;
 fnu = (2.0e0 - fid) / 3.0e0;
 nz = (0, _zbinu.zbinu)(ztar, ztai, fnu, kode, 2, cyr, cyi, rl, fnul, tol, elim, alim);
 cyr[0] = cyr[0] * z3r;
 cyi[0] = cyi[0] * z3r;
 cyr[1] = cyr[1] * z3r;
 cyi[1] = cyi[1] * z3r;
 // -----------------------------------------------------------------------
 // backward recur one step for orders -Math.trunc(1/3) or -Math.trunc(2/3)
 // -----------------------------------------------------------------------

 var _zdiv = (0, _zdiv3.zdiv)(cyr[0], cyi[0], ztar, ztai);

 var _zdiv2 = _slicedToArray(_zdiv, 2);

 str = _zdiv2[0];
 sti = _zdiv2[1];

 s2r = (fnu + fnu) * str + cyr[1];
 s2i = (fnu + fnu) * sti + cyi[1];
 aa = fmr * (fnu - 1.0e0);
 str = Math.cos(aa);
 sti = Math.sin(aa);
 s1r = coef * (s1r + s2r * str - s2i * sti);
 s1i = coef * (s1i + s2r * sti + s2i * str);
 if (id === 1) {
 goToLabel = 120;break;
 }
 str = csqr * s1r - csqi * s1i;
 s1i = csqr * s1i + csqi * s1r;
 s1r = str;
 bir = s1r / sfac;
 bii = s1i / sfac;
 break mainExecutionLoop;
 case 120:

 str = zr * s1r - zi * s1i;
 s1i = zr * s1i + zi * s1r;
 s1r = str;
 bir = s1r / sfac;
 bii = s1i / sfac;
 break mainExecutionLoop;
 case 130:
 aa = c1 * (1.0e0 - fid) + fid * c2;
 bir = aa;
 bii = 0.0e0;
 break mainExecutionLoop;
 case 190:

 ierr = 2;
 nz = 0;
 break mainExecutionLoop;
 case 200:

 if (nz === -1) {
 goToLabel = 190;break;
 }
 nz = 0;
 ierr = 5;
 break mainExecutionLoop;
 case 260:

 ierr = 4;
 nz = 0;

 default:
 break mainExecutionLoop;
 }
 }
 return [bir, bii, nz, ierr];
}
},{"../../utils/fortran-utils/d1mach.js":91,"../../utils/fortran-utils/i1mach.js":92,"./zabs.js":11,"./zbinu.js":21,"./zdiv.js":26,"./zsqrt.js":36}],23:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});

var _slicedToArray = function () { function sliceIterator(arr, i) { var _arr = []; var _n = true; var _d = false; var _e = undefined; try { for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) { _arr.push(_s.value); if (i && _arr.length === i) break; } } catch (err) { _d = true; _e = err; } finally { try { if (!_n && _i["return"]) _i["return"](); } finally { if (_d) throw _e; } } return _arr; } return function (arr, i) { if (Array.isArray(arr)) { return arr; } else if (Symbol.iterator in Object(arr)) { return sliceIterator(arr, i); } else { throw new TypeError("Invalid attempt to destructure non-iterable instance"); } }; }(); /* eslint-disable no-labels */
/* eslint-disable no-fallthrough */
// SUBROUTINE ZBKNU(ZR, ZI, FNU, KODE, N, YR, YI, NZ, TOL, ELIM, ALIM)
// ***BEGIN PROLOGUE ZBKNU
// ***REFER TO ZBESI,ZBESK,ZAIRY,ZBESH
//
// ZBKNU COMPUTES THE K BESSEL FUNCTION IN THE RIGHT HALF Z PLANE.
//
// ***ROUTINES CALLED DGAMLN,I1MACH,D1MACH,ZKSCL,ZSHCH,ZUCHK,AZABS,ZDIV,
// AZEXP,AZLOG,ZMLT,AZSQRT
// ***END PROLOGUE ZBKNU
//


exports.zbknu = zbknu;

var _dgamln = require('./dgamln.js');

var _i1mach = require('../../utils/fortran-utils/i1mach.js');

var _d1mach = require('../../utils/fortran-utils/d1mach.js');

var _zkscl = require('./zkscl.js');

var _zshch3 = require('./zshch.js');

var _zuchk = require('./zuchk.js');

var _zabs = require('./zabs.js');

var _zdiv7 = require('./zdiv.js');

var _zexp = require('./zexp.js');

var _zlog = require('./zlog.js');

var _zmlt23 = require('./zmlt.js');

var _zsqrt = require('./zsqrt.js');

function zbknu(zr, zi, fnu, kode, n, yr, yi, tol, elim, alim) {
 var aa = void 0,
 ak = void 0,
 ascle = void 0,
 a1 = void 0,
 a2 = void 0,
 bb = void 0,
 bk = void 0,
 bry = void 0,
 caz = void 0,
 cbi = void 0,
 cbr = void 0,
 cc = void 0,
 cchi = void 0,
 cchr = void 0,
 cki = void 0,
 ckr = void 0,
 coefi = void 0,
 coefr = void 0,
 conei = void 0,
 coner = void 0,
 crscr = void 0,
 csclr = void 0,
 cshi = void 0,
 cshr = void 0,
 csi = void 0,
 csr = void 0,
 csrr = void 0,
 cssr = void 0,
 ctwor = void 0,
 czeroi = void 0,
 czeror = void 0,
 czi = void 0,
 czr = void 0,
 dnu = void 0,
 dnu2 = void 0,
 dpi = void 0,
 etest = void 0,
 fc = void 0,
 fhs = void 0,
 fi = void 0,
 fk = void 0,
 fks = void 0,
 fmui = void 0,
 fmur = void 0,
 fpi = void 0,
 fr = void 0,
 g1 = void 0,
 g2 = void 0,
 hpi = void 0,
 pi = void 0,
 pr = void 0,
 pti = void 0,
 ptr = void 0,
 p1i = void 0,
 p1r = void 0,
 p2i = void 0,
 p2m = void 0,
 p2r = void 0,
 qi = void 0,
 qr = void 0,
 rak = void 0,
 rcaz = void 0,
 rthpi = void 0,
 rzi = void 0,
 rzr = void 0,
 r1 = void 0,
 s = void 0,
 smui = void 0,
 smur = void 0,
 spi = void 0,
 sti = void 0,
 str = void 0,
 s1i = void 0,
 s1r = void 0,
 s2i = void 0,
 s2r = void 0,
 tm = void 0,
 tth = void 0,
 t1 = void 0,
 t2 = void 0,
 elm = void 0,
 celmr = void 0,
 zdr = void 0,
 zdi = void 0,
 as = void 0,
 alas = void 0,
 helim = void 0,
 cyr = void 0,
 cyi = void 0,
 i = void 0,
 iflag = void 0,
 inu = void 0,
 k = void 0,
 kflag = void 0,
 kk = void 0,
 kmax = void 0,
 koded = void 0,
 nz = void 0,
 idum = void 0,
 j = void 0,
 ic = void 0,
 inub = void 0,
 nw = void 0;

 cssr = new Array(3);
 csrr = new Array(3);
 bry = new Array(3);
 cyr = new Array(2);
 cyi = new Array(2);

 kmax = 30;

 czeror = 0.0;
 czeroi = 0.0;
 coner = 1.0;
 conei = 0.0;
 ctwor = 2.0;
 r1 = 2.0;
 dpi = 3.14159265358979324;
 rthpi = 1.25331413731550025;
 spi = 1.90985931710274403;
 hpi = 1.57079632679489662;
 fpi = 1.89769999331517738;
 tth = 6.66666666666666666e-01;


 cc = [5.77215664901532861e-01, -4.20026350340952355e-02, -4.21977345555443367e-02, 7.21894324666309954e-03, -2.15241674114950973e-04, -2.01348547807882387e-05, 1.13302723198169588e-06, 6.11609510448141582e-09];

 var goToLabel = 0;
 mainExecutionLoop: while (true) {
 switch (goToLabel) {
 case 0:
 caz = (0, _zabs.azabs)(zr, zi);
 csclr = 1.0 / tol;
 crscr = tol;
 cssr[0] = csclr;
 cssr[1] = 1.0;
 cssr[2] = crscr;
 csrr[0] = crscr;
 csrr[1] = 1.0;
 csrr[2] = csclr;
 bry[0] = 1.0e+3 * (0, _d1mach.d1mach)(1) / tol;
 bry[1] = 1.0 / bry[0];
 bry[2] = (0, _d1mach.d1mach)(2);
 nz = 0;
 iflag = 0;
 koded = kode;
 rcaz = 1.0 / caz;
 str = zr * rcaz;
 sti = -zi * rcaz;
 rzr = (str + str) * rcaz;
 rzi = (sti + sti) * rcaz;
 inu = Math.trunc(fnu + 0.5);
 dnu = fnu - inu;
 if (Math.abs(dnu) === 0.5) {
 goToLabel = 110;break;
 }
 dnu2 = 0.0;
 if (Math.abs(dnu) > tol) dnu2 = dnu * dnu;
 if (caz > r1) {
 goToLabel = 110;break;
 }
 // c-----------------------------------------------------------------------
 // c series for cabs(z) <= r1
 // c-----------------------------------------------------------------------
 fc = 1.0;

 var _azlog = (0, _zlog.azlog)(rzr, rzi);

 var _azlog2 = _slicedToArray(_azlog, 3);

 smur = _azlog2[0];
 smui = _azlog2[1];
 idum = _azlog2[2];

 fmur = smur * dnu;
 fmui = smui * dnu;

 var _zshch = (0, _zshch3.zshch)(fmur, fmui);

 var _zshch2 = _slicedToArray(_zshch, 4);

 cshr = _zshch2[0];
 cshi = _zshch2[1];
 cchr = _zshch2[2];
 cchi = _zshch2[3];

 if (dnu === 0.0) {
 goToLabel = 10;break;
 }
 fc = dnu * dpi;
 fc = fc / Math.sin(fc);
 smur = cshr / dnu;
 smui = cshi / dnu;
 case 10:
 a2 = 1.0 + dnu;
 // c-----------------------------------------------------------------------
 // c gam(1-z)*gam(1+z)=pi*z/sin(pi*z), t1=1/gam(1-dnu), t2=1/gam(1+dnu)
 // c-----------------------------------------------------------------------
 t2 = Math.exp(-(0, _dgamln.dgamln)(a2, idum));
 t1 = 1.0 / (t2 * fc);
 if (Math.abs(dnu) > 0.1) {
 goToLabel = 40;break;
 }
 // c-----------------------------------------------------------------------
 // c series for f0 to resolve indeterminacy for small abs(dnu)
 // c-----------------------------------------------------------------------
 ak = 1.0;
 s = cc[0];
 // do 20 k=2,8
 for (k = 2; k <= 8; k++) {
 ak = ak * dnu2;
 tm = cc[k - 1] * ak;
 s = s + tm;
 if (Math.abs(tm) < tol) break;
 }
 // 20 continue
 g1 = -s;
 goToLabel = 50;break;
 case 40:
 g1 = (t1 - t2) / (dnu + dnu);
 case 50:
 g2 = (t1 + t2) * 0.5;
 fr = fc * (cchr * g1 + smur * g2);
 fi = fc * (cchi * g1 + smui * g2);

 var _azexp = (0, _zexp.azexp)(fmur, fmui);

 var _azexp2 = _slicedToArray(_azexp, 2);

 str = _azexp2[0];
 sti = _azexp2[1];

 pr = 0.5 * str / t2;
 pi = 0.5 * sti / t2;

 var _zdiv = (0, _zdiv7.zdiv)(0.5, 0.0, str, sti);

 var _zdiv2 = _slicedToArray(_zdiv, 2);

 ptr = _zdiv2[0];
 pti = _zdiv2[1];

 qr = ptr / t1;
 qi = pti / t1;
 s1r = fr;
 s1i = fi;
 s2r = pr;
 s2i = pi;
 ak = 1.0;
 a1 = 1.0;
 ckr = coner;
 cki = conei;
 bk = 1.0 - dnu2;
 if (inu > 0 || n > 1) {
 goToLabel = 80;break;
 }
 // c-----------------------------------------------------------------------
 // c generate k(fnu,z), 0.0 <= fnu < 0.5 and n=1
 // c-----------------------------------------------------------------------
 if (caz < tol) {
 goToLabel = 70;break;
 }

 var _zmlt = (0, _zmlt23.zmlt)(zr, zi, zr, zi);

 var _zmlt2 = _slicedToArray(_zmlt, 2);

 czr = _zmlt2[0];
 czi = _zmlt2[1];

 czr = 0.25 * czr;
 czi = 0.25 * czi;
 t1 = 0.25 * caz * caz;
 case 60:
 fr = (fr * ak + pr + qr) / bk;
 fi = (fi * ak + pi + qi) / bk;
 str = 1.0 / (ak - dnu);
 pr = pr * str;
 pi = pi * str;
 str = 1.0 / (ak + dnu);
 qr = qr * str;
 qi = qi * str;
 str = ckr * czr - cki * czi;
 rak = 1.0 / ak;
 cki = (ckr * czi + cki * czr) * rak;
 ckr = str * rak;
 s1r = ckr * fr - cki * fi + s1r;
 s1i = ckr * fi + cki * fr + s1i;
 a1 = a1 * t1 * rak;
 bk = bk + ak + ak + 1.0;
 ak = ak + 1.0;
 if (a1 > tol) {
 goToLabel = 60;break;
 }
 case 70:
 yr[0] = s1r;
 yi[0] = s1i;
 if (koded === 1) break mainExecutionLoop;

 var _azexp3 = (0, _zexp.azexp)(zr, zi);

 var _azexp4 = _slicedToArray(_azexp3, 2);

 str = _azexp4[0];
 sti = _azexp4[1];

 var _zmlt3 = (0, _zmlt23.zmlt)(s1r, s1i, str, sti);

 var _zmlt4 = _slicedToArray(_zmlt3, 2);

 yr[0] = _zmlt4[0];
 yi[0] = _zmlt4[1];

 break mainExecutionLoop;
 // c-----------------------------------------------------------------------
 // c generate k(dnu,z) and k(dnu+1,z) for forward recurrence
 // c-----------------------------------------------------------------------
 case 80:
 if (caz < tol) {
 goToLabel = 100;break;
 }

 var _zmlt5 = (0, _zmlt23.zmlt)(zr, zi, zr, zi);

 var _zmlt6 = _slicedToArray(_zmlt5, 2);

 czr = _zmlt6[0];
 czi = _zmlt6[1];

 czr = 0.25 * czr;
 czi = 0.25 * czi;
 t1 = 0.25 * caz * caz;
 case 90:
 fr = (fr * ak + pr + qr) / bk;
 fi = (fi * ak + pi + qi) / bk;
 str = 1.0 / (ak - dnu);
 pr = pr * str;
 pi = pi * str;
 str = 1.0 / (ak + dnu);
 qr = qr * str;
 qi = qi * str;
 str = ckr * czr - cki * czi;
 rak = 1.0 / ak;
 cki = (ckr * czi + cki * czr) * rak;
 ckr = str * rak;
 s1r = ckr * fr - cki * fi + s1r;
 s1i = ckr * fi + cki * fr + s1i;
 str = pr - fr * ak;
 sti = pi - fi * ak;
 s2r = ckr * str - cki * sti + s2r;
 s2i = ckr * sti + cki * str + s2i;
 a1 = a1 * t1 * rak;
 bk = bk + ak + ak + 1.0;
 ak = ak + 1.0;
 if (a1 > tol) {
 goToLabel = 90;break;
 }
 case 100:
 kflag = 2;
 a1 = fnu + 1.0;
 ak = a1 * Math.abs(smur);
 if (ak > alim) kflag = 3;
 str = cssr[kflag - 1];
 p2r = s2r * str;
 p2i = s2i * str;

 var _zmlt7 = (0, _zmlt23.zmlt)(p2r, p2i, rzr, rzi);

 var _zmlt8 = _slicedToArray(_zmlt7, 2);

 s2r = _zmlt8[0];
 s2i = _zmlt8[1];

 s1r = s1r * str;
 s1i = s1i * str;
 if (koded === 1) {
 goToLabel = 210;break;
 }

 var _azexp5 = (0, _zexp.azexp)(zr, zi);

 var _azexp6 = _slicedToArray(_azexp5, 2);

 fr = _azexp6[0];
 fi = _azexp6[1];

 var _zmlt9 = (0, _zmlt23.zmlt)(s1r, s1i, fr, fi);

 var _zmlt10 = _slicedToArray(_zmlt9, 2);

 s1r = _zmlt10[0];
 s1i = _zmlt10[1];

 var _zmlt11 = (0, _zmlt23.zmlt)(s2r, s2i, fr, fi);

 var _zmlt12 = _slicedToArray(_zmlt11, 2);

 s2r = _zmlt12[0];
 s2i = _zmlt12[1];

 goToLabel = 210;break;
 // c-----------------------------------------------------------------------
 // c iflag=0 means no underflow occurred
 // c iflag=1 means an underflow occurred- computation proceeds with
 // c koded=2 and a test for on scale values is made during forward
 // c recursion
 // c-----------------------------------------------------------------------
 case 110:
 var _azsqrt = (0, _zsqrt.azsqrt)(zr, zi);

 var _azsqrt2 = _slicedToArray(_azsqrt, 2);

 str = _azsqrt2[0];
 sti = _azsqrt2[1];

 var _zdiv3 = (0, _zdiv7.zdiv)(rthpi, czeroi, str, sti);

 var _zdiv4 = _slicedToArray(_zdiv3, 2);

 coefr = _zdiv4[0];
 coefi = _zdiv4[1];

 kflag = 2;
 if (koded === 2) {
 goToLabel = 120;break;
 }
 if (zr > alim) {
 goToLabel = 290;break;
 }
 // c blank line
 str = Math.exp(-zr) * cssr[kflag - 1];
 sti = -str * Math.sin(zi);
 str = str * Math.cos(zi);

 var _zmlt13 = (0, _zmlt23.zmlt)(coefr, coefi, str, sti);

 var _zmlt14 = _slicedToArray(_zmlt13, 2);

 coefr = _zmlt14[0];
 coefi = _zmlt14[1];

 case 120:
 if (Math.abs(dnu) === 0.5) {
 goToLabel = 300;break;
 }
 // c-----------------------------------------------------------------------
 // c miller algorithm for cabs(z) > r1
 // c-----------------------------------------------------------------------
 ak = Math.cos(dpi * dnu);
 ak = Math.abs(ak);
 if (ak === czeror) {
 goToLabel = 300;break;
 }
 fhs = Math.abs(0.25 - dnu2);
 if (fhs === czeror) {
 goToLabel = 300;break;
 }
 // c-----------------------------------------------------------------------
 // c compute r2=f(e). if cabs(z) >= r2, use forward recurrence to
 // c determine the backward index k. r2=f(e) is a straight line on
 // c 12 <= e <= 60. e is computed from 2**(-e)=b**(1-i1mach(14))=
 // c tol where b is the base of the arithmetic.
 // c-----------------------------------------------------------------------
 t1 = (0, _i1mach.i1mach)(14) - 1;
 t1 = t1 * (0, _d1mach.d1mach)(5) * 3.321928094;
 t1 = Math.max(t1, 12.0);
 t1 = Math.min(t1, 60.0);
 t2 = tth * t1 - 6.0;
 if (zr !== 0.0) {
 goToLabel = 130;break;
 }
 t1 = hpi;
 goToLabel = 140;break;
 case 130:
 t1 = Math.atan(zi / zr);
 t1 = Math.abs(t1);
 case 140:
 if (t2 > caz) {
 goToLabel = 170;break;
 }
 // c-----------------------------------------------------------------------
 // c forward recurrence loop when cabs(z) >= r2
 // c-----------------------------------------------------------------------
 etest = ak / (dpi * caz * tol);
 fk = coner;
 if (etest < coner) {
 goToLabel = 180;break;
 }
 fks = ctwor;
 ckr = caz + caz + ctwor;
 p1r = czeror;
 p2r = coner;
 // do 150 i=1,kmax
 for (i = 1; i <= kmax; i++) {
 ak = fhs / fks;
 cbr = ckr / (fk + coner);
 ptr = p2r;
 p2r = cbr * p2r - p1r * ak;
 p1r = ptr;
 ckr = ckr + ctwor;
 fks = fks + fk + fk + ctwor;
 fhs = fhs + fk + fk;
 fk = fk + coner;
 str = Math.abs(p2r) * fk;
 if (etest < str) {
 goToLabel = 160;break;
 }
 }
 // 150 continue
 if (goToLabel === 160) break;
 goToLabel = 310;break;
 case 160:
 fk = fk + spi * t1 * Math.sqrt(t2 / caz);
 fhs = Math.abs(0.25 - dnu2);
 goToLabel = 180;break;
 case 170:
 // c-----------------------------------------------------------------------
 // c compute backward index k for cabs(z) < r2
 // c-----------------------------------------------------------------------
 a2 = Math.sqrt(caz);
 ak = fpi * ak / (tol * Math.sqrt(a2));
 aa = 3.0 * t1 / (1.0 + caz);
 bb = 14.7 * t1 / (28.0 + caz);
 ak = (Math.log(ak) + caz * Math.cos(aa) / (1.0 + 0.008 * caz)) / Math.cos(bb);
 fk = 0.12125 * ak * ak / caz + 1.5;
 case 180:
 // c-----------------------------------------------------------------------
 // c backward recurrence loop for miller algorithm
 // c-----------------------------------------------------------------------
 k = Math.trunc(fk);
 fk = k;
 fks = fk * fk;
 p1r = czeror;
 p1i = czeroi;
 p2r = tol;
 p2i = czeroi;
 csr = p2r;
 csi = p2i;
 // do 190 i=1,k
 for (i = 1; i <= k; i++) {
 a1 = fks - fk;
 ak = (fks + fk) / (a1 + fhs);
 rak = 2.0 / (fk + coner);
 cbr = (fk + zr) * rak;
 cbi = zi * rak;
 ptr = p2r;
 pti = p2i;
 p2r = (ptr * cbr - pti * cbi - p1r) * ak;
 p2i = (pti * cbr + ptr * cbi - p1i) * ak;
 p1r = ptr;
 p1i = pti;
 csr = csr + p2r;
 csi = csi + p2i;
 fks = a1 - fk + coner;
 fk = fk - coner;
 }
 // 190 continue
 // c-----------------------------------------------------------------------
 // c compute (p2/cs)=(p2/cabs(cs))*(conjg(cs)/cabs(cs)) for better
 // c scaling
 // c-----------------------------------------------------------------------
 tm = (0, _zabs.azabs)(csr, csi);
 ptr = 1.0 / tm;
 s1r = p2r * ptr;
 s1i = p2i * ptr;
 csr = csr * ptr;
 csi = -csi * ptr;

 var _zmlt15 = (0, _zmlt23.zmlt)(coefr, coefi, s1r, s1i);

 var _zmlt16 = _slicedToArray(_zmlt15, 2);

 str = _zmlt16[0];
 sti = _zmlt16[1];

 var _zmlt17 = (0, _zmlt23.zmlt)(str, sti, csr, csi);

 var _zmlt18 = _slicedToArray(_zmlt17, 2);

 s1r = _zmlt18[0];
 s1i = _zmlt18[1];

 if (inu > 0 || n > 1) {
 goToLabel = 200;break;
 }
 zdr = zr;
 zdi = zi;
 if (iflag === 1) {
 goToLabel = 270;break;
 }
 goToLabel = 240;break;
 case 200:
 // c-----------------------------------------------------------------------
 // c compute p1/p2=(p1/cabs(p2)*conjg(p2)/cabs(p2) for scaling
 // c-----------------------------------------------------------------------
 tm = (0, _zabs.azabs)(p2r, p2i);
 ptr = 1.0 / tm;
 p1r = p1r * ptr;
 p1i = p1i * ptr;
 p2r = p2r * ptr;
 p2i = -p2i * ptr;

 var _zmlt19 = (0, _zmlt23.zmlt)(p1r, p1i, p2r, p2i);

 var _zmlt20 = _slicedToArray(_zmlt19, 2);

 ptr = _zmlt20[0];
 pti = _zmlt20[1];

 str = dnu + 0.5 - ptr;
 sti = -pti;

 var _zdiv5 = (0, _zdiv7.zdiv)(str, sti, zr, zi);

 var _zdiv6 = _slicedToArray(_zdiv5, 2);

 str = _zdiv6[0];
 sti = _zdiv6[1];

 str = str + 1.0;

 var _zmlt21 = (0, _zmlt23.zmlt)(str, sti, s1r, s1i);

 var _zmlt22 = _slicedToArray(_zmlt21, 2);

 s2r = _zmlt22[0];
 s2i = _zmlt22[1];

 // c-----------------------------------------------------------------------
 // c forward recursion on the three term recursion with relation with
 // c scaling near exponent extremes on kflag=1 or kflag=3
 // c-----------------------------------------------------------------------
 case 210:
 str = dnu + 1.0;
 ckr = str * rzr;
 cki = str * rzi;
 if (n === 1) inu = inu - 1;
 if (inu > 0) {
 goToLabel = 220;break;
 }
 if (n > 1) {
 goToLabel = 215;break;
 }
 s1r = s2r;
 s1i = s2i;
 case 215:
 zdr = zr;
 zdi = zi;
 if (iflag === 1) {
 goToLabel = 270;break;
 }
 goToLabel = 240;break;
 case 220:
 inub = 1;
 if (iflag === 1) {
 goToLabel = 261;break;
 }
 case 225:
 p1r = csrr[kflag - 1];
 ascle = bry[kflag - 1];
 // do 230 i=inub,inu
 for (i = inub; i <= inu; i++) {
 str = s2r;
 sti = s2i;
 s2r = ckr * str - cki * sti + s1r;
 s2i = ckr * sti + cki * str + s1i;
 s1r = str;
 s1i = sti;
 ckr = ckr + rzr;
 cki = cki + rzi;
 if (kflag >= 3) continue;
 p2r = s2r * p1r;
 p2i = s2i * p1r;
 str = Math.abs(p2r);
 sti = Math.abs(p2i);
 p2m = Math.max(str, sti);
 if (p2m <= ascle) continue;
 kflag = kflag + 1;
 ascle = bry[kflag - 1];
 s1r = s1r * p1r;
 s1i = s1i * p1r;
 s2r = p2r;
 s2i = p2i;
 str = cssr[kflag - 1];
 s1r = s1r * str;
 s1i = s1i * str;
 s2r = s2r * str;
 s2i = s2i * str;
 p1r = csrr[kflag - 1];
 }
 // 230 continue
 if (n !== 1) {
 goToLabel = 240;break;
 }
 s1r = s2r;
 s1i = s2i;
 case 240:
 str = csrr[kflag - 1];
 yr[0] = s1r * str;
 yi[0] = s1i * str;
 if (n === 1) break mainExecutionLoop;
 yr[1] = s2r * str;
 yi[1] = s2i * str;
 if (n === 2) break mainExecutionLoop;
 kk = 2;
 case 250:
 kk = kk + 1;
 if (kk > n) break mainExecutionLoop;
 p1r = csrr[kflag - 1];
 ascle = bry[kflag - 1];
 // do 260 i=kk,n
 for (i = kk; i <= n; i++) {
 p2r = s2r;
 p2i = s2i;
 s2r = ckr * p2r - cki * p2i + s1r;
 s2i = cki * p2r + ckr * p2i + s1i;
 s1r = p2r;
 s1i = p2i;
 ckr = ckr + rzr;
 cki = cki + rzi;
 p2r = s2r * p1r;
 p2i = s2i * p1r;
 yr[i - 1] = p2r;
 yi[i - 1] = p2i;
 if (kflag >= 3) continue;
 str = Math.abs(p2r);
 sti = Math.abs(p2i);
 p2m = Math.max(str, sti);
 if (p2m <= ascle) continue;
 kflag = kflag + 1;
 ascle = bry[kflag - 1];
 s1r = s1r * p1r;
 s1i = s1i * p1r;
 s2r = p2r;
 s2i = p2i;
 str = cssr[kflag - 1];
 s1r = s1r * str;
 s1i = s1i * str;
 s2r = s2r * str;
 s2i = s2i * str;
 p1r = csrr[kflag - 1];
 }
 // 260 continue
 break mainExecutionLoop;
 // c-----------------------------------------------------------------------
 // c iflag=1 cases, forward recurrence on scaled values on underflow
 // c-----------------------------------------------------------------------
 case 261:
 helim = 0.5 * elim;
 elm = Math.exp(-elim);
 celmr = elm;
 ascle = bry[0];
 zdr = zr;
 zdi = zi;
 ic = -1;
 j = 2;
 // do 262 i=1,inu
 for (i = 1; i <= inu; i++) {
 str = s2r;
 sti = s2i;
 s2r = str * ckr - sti * cki + s1r;
 s2i = sti * ckr + str * cki + s1i;
 s1r = str;
 s1i = sti;
 ckr = ckr + rzr;
 cki = cki + rzi;
 as = (0, _zabs.azabs)(s2r, s2i);
 alas = Math.log(as);
 p2r = -zdr + alas;
 if (p2r < -elim) {
 // go to 263
 } else {
 var _azlog3 = (0, _zlog.azlog)(s2r, s2i);

 var _azlog4 = _slicedToArray(_azlog3, 3);

 str = _azlog4[0];
 sti = _azlog4[1];
 idum = _azlog4[2];

 p2r = -zdr + str;
 p2i = -zdi + sti;
 p2m = Math.exp(p2r) / tol;
 p1r = p2m * Math.cos(p2i);
 p1i = p2m * Math.sin(p2i);
 nw = (0, _zuchk.zuchk)(p1r, p1i, ascle, tol);
 if (nw !== 0) {
 // go to 263
 } else {
 j = 3 - j;
 cyr[j - 1] = p1r;
 cyi[j - 1] = p1i;
 if (ic === i - 1) {
 goToLabel = 264;break;
 }
 ic = i;
 continue;
 }
 }
 // 263 continue
 if (alas < helim) continue;
 zdr = zdr - elim;
 s1r = s1r * celmr;
 s1i = s1i * celmr;
 s2r = s2r * celmr;
 s2i = s2i * celmr;
 }
 // 262 continue
 if (goToLabel === 264) break;
 if (n !== 1) {
 goToLabel = 270;break;
 }
 s1r = s2r;
 s1i = s2i;
 goToLabel = 270;break;
 case 264:
 kflag = 1;
 inub = i + 1;
 s2r = cyr[j - 1];
 s2i = cyi[j - 1];
 j = 3 - j;
 s1r = cyr[j - 1];
 s1i = cyi[j - 1];
 if (inub <= inu) {
 goToLabel = 225;break;
 }
 if (n !== 1) {
 goToLabel = 240;break;
 }
 s1r = s2r;
 s1i = s2i;
 goToLabel = 240;break;
 case 270:
 yr[0] = s1r;
 yi[0] = s1i;
 if (n === 1) {
 goToLabel = 280;break;
 }
 yr[1] = s2r;
 yi[1] = s2i;
 case 280:
 ascle = bry[0];
 nz = (0, _zkscl.zkscl)(zdr, zdi, fnu, n, yr, yi, rzr, rzi, ascle, tol, elim);
 inu = n - nz;
 if (inu <= 0) break mainExecutionLoop;
 kk = nz + 1;
 s1r = yr[kk - 1];
 s1i = yi[kk - 1];
 yr[kk - 1] = s1r * csrr[0];
 yi[kk - 1] = s1i * csrr[0];
 if (inu === 1) break mainExecutionLoop;
 kk = nz + 2;
 s2r = yr[kk - 1];
 s2i = yi[kk - 1];
 yr[kk - 1] = s2r * csrr[0];
 yi[kk - 1] = s2i * csrr[0];
 if (inu === 2) break mainExecutionLoop;
 t2 = fnu + (kk - 1);
 ckr = t2 * rzr;
 cki = t2 * rzi;
 kflag = 1;
 goToLabel = 250;break;
 case 290:
 // c-----------------------------------------------------------------------
 // c scale by Math.exp(z), iflag = 1 cases
 // c-----------------------------------------------------------------------
 koded = 2;
 iflag = 1;
 kflag = 2;
 goToLabel = 120;break;
 // c-----------------------------------------------------------------------
 // c fnu=half odd integer case, dnu=-0.5
 // c-----------------------------------------------------------------------
 case 300:
 s1r = coefr;
 s1i = coefi;
 s2r = coefr;
 s2i = coefi;
 goToLabel = 210;break;
 // c
 // c
 case 310:
 nz = -2;
 default:
 break mainExecutionLoop;
 }
 }

 return nz;
}
},{"../../utils/fortran-utils/d1mach.js":91,"../../utils/fortran-utils/i1mach.js":92,"./dgamln.js":10,"./zabs.js":11,"./zdiv.js":26,"./zexp.js":27,"./zkscl.js":28,"./zlog.js":29,"./zmlt.js":31,"./zshch.js":35,"./zsqrt.js":36,"./zuchk.js":37}],24:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});

var _slicedToArray = function () { function sliceIterator(arr, i) { var _arr = []; var _n = true; var _d = false; var _e = undefined; try { for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) { _arr.push(_s.value); if (i && _arr.length === i) break; } } catch (err) { _d = true; _e = err; } finally { try { if (!_n && _i["return"]) _i["return"](); } finally { if (_d) throw _e; } } return _arr; } return function (arr, i) { if (Array.isArray(arr)) { return arr; } else if (Symbol.iterator in Object(arr)) { return sliceIterator(arr, i); } else { throw new TypeError("Invalid attempt to destructure non-iterable instance"); } }; }(); /* eslint-disable no-labels */
/* eslint-disable no-fallthrough */
// SUBROUTINE ZBUNI(ZR, ZI, FNU, KODE, N, YR, YI, NZ, NUI, NLAST,
// * FNUL, TOL, ELIM, ALIM)
// ***BEGIN PROLOGUE ZBUNI
// ***REFER TO ZBESI,ZBESK
//
// ZBUNI COMPUTES THE I BESSEL FUNCTION FOR LARGE CABS(Z).GT.
// FNUL AND FNU+N-1.LT.FNUL. THE ORDER IS INCREASED FROM
// FNU+N-1 GREATER THAN FNUL BY ADDING NUI AND COMPUTING
// ACCORDING TO THE UNIFORM ASYMPTOTIC EXPANSION FOR I(FNU,Z)
// ON IFORM=1 AND THE EXPANSION FOR J(FNU,Z) ON IFORM=2
//
// ***ROUTINES CALLED ZUNI1,ZUNI2,AZABS,D1MACH
// ***END PROLOGUE ZBUNI


exports.zbuni = zbuni;

var _zuni9 = require('./zuni1.js');

var _zuni10 = require('./zuni2.js');

var _zabs = require('./zabs.js');

var _d1mach = require('../../utils/fortran-utils/d1mach.js');

function zbuni(zr, zi, fnu, kode, n, yr, yi, nui, nlast, fnul, tol, elim, alim) {
 var ax = void 0,
 ay = void 0,
 csclr = void 0,
 cscrr = void 0,
 cyi = void 0,
 cyr = void 0,
 dfnu = void 0,
 fnui = void 0,
 gnu = void 0,
 raz = void 0,
 rzi = void 0,
 rzr = void 0,
 sti = void 0,
 str = void 0,
 s1i = void 0,
 s1r = void 0,
 s2i = void 0,
 s2r = void 0,
 ascle = void 0,
 bry = void 0,
 c1r = void 0,
 c1i = void 0,
 c1m = void 0,
 i = void 0,
 iflag = void 0,
 iform = void 0,
 k = void 0,
 nl = void 0,
 nw = void 0,
 nz = void 0;

 cyr = new Array(2);
 cyi = new Array(2);
 bry = new Array(3);

 var goToLabel = 0;
 mainExecutionLoop: while (true) {
 switch (goToLabel) {
 case 0:
 nz = 0;
 ax = Math.abs(zr) * 1.7321;
 ay = Math.abs(zi);
 iform = 1;
 if (ay > ax) iform = 2;
 if (nui === 0) {
 goToLabel = 60;break;
 }
 fnui = nui;
 dfnu = fnu + (n - 1);
 gnu = dfnu + fnui;
 if (iform === 2) {
 goToLabel = 10;break;
 }
 // c-----------------------------------------------------------------------
 // c asymptotic expansion for i(fnu,z) for large fnu applied in
 // c -pi/3 <= arg(z) <= pi/3
 // c-----------------------------------------------------------------------

 var _zuni = (0, _zuni9.zuni1)(zr, zi, gnu, kode, 2, cyr, cyi, fnul, tol, elim, alim);

 var _zuni2 = _slicedToArray(_zuni, 2);

 nw = _zuni2[0];
 nlast = _zuni2[1];

 goToLabel = 20;break;
 case 10:
 var _zuni3 = (0, _zuni10.zuni2)(zr, zi, gnu, kode, 2, cyr, cyi, fnul, tol, elim, alim);
 // c-----------------------------------------------------------------------
 // c asymptotic expansion for j(fnu,z*exp(m*hpi)) for large fnu
 // c applied in pi/3 < abs(arg(z)) <= pi/2 where m=+i or -i
 // c and hpi=pi/2
 // c-----------------------------------------------------------------------


 var _zuni4 = _slicedToArray(_zuni3, 2);

 nw = _zuni4[0];
 nlast = _zuni4[1];

 case 20:
 if (nw < 0) {
 goToLabel = 50;break;
 }
 if (nw !== 0) {
 goToLabel = 90;break;
 }
 str = (0, _zabs.azabs)(cyr[0], cyi[0]);
 // c----------------------------------------------------------------------
 // c scale backward recurrence, bry(3) is defined but never used
 // c----------------------------------------------------------------------
 bry[0] = 1.0e+3 * (0, _d1mach.d1mach)(1) / tol;
 bry[1] = 1.0 / bry[0];
 bry[2] = bry[1];
 iflag = 2;
 ascle = bry[1];
 csclr = 1.0;
 if (str > bry[0]) {
 goToLabel = 21;break;
 }
 iflag = 1;
 ascle = bry[0];
 csclr = 1.0 / tol;
 goToLabel = 25;break;
 case 21:
 if (str < bry[1]) {
 goToLabel = 25;break;
 }
 iflag = 3;
 ascle = bry[2];
 csclr = tol;
 case 25:
 cscrr = 1.0 / csclr;
 s1r = cyr[1] * csclr;
 s1i = cyi[1] * csclr;
 s2r = cyr[0] * csclr;
 s2i = cyi[0] * csclr;
 raz = 1.0 / (0, _zabs.azabs)(zr, zi);
 str = zr * raz;
 sti = -zi * raz;
 rzr = (str + str) * raz;
 rzi = (sti + sti) * raz;
 // do 30 i=1,nui
 for (i = 1; i <= nui; i++) {
 str = s2r;
 sti = s2i;
 s2r = (dfnu + fnui) * (rzr * str - rzi * sti) + s1r;
 s2i = (dfnu + fnui) * (rzr * sti + rzi * str) + s1i;
 s1r = str;
 s1i = sti;
 fnui = fnui - 1.0;
 if (iflag >= 3) continue;
 str = s2r * cscrr;
 sti = s2i * cscrr;
 c1r = Math.abs(str);
 c1i = Math.abs(sti);
 c1m = Math.max(c1r, c1i);
 if (c1m <= ascle) continue;
 iflag = iflag + 1;
 ascle = bry[iflag - 1];
 s1r = s1r * cscrr;
 s1i = s1i * cscrr;
 s2r = str;
 s2i = sti;
 csclr = csclr * tol;
 cscrr = 1.0 / csclr;
 s1r = s1r * csclr;
 s1i = s1i * csclr;
 s2r = s2r * csclr;
 s2i = s2i * csclr;
 }
 // 30 continue
 yr[n - 1] = s2r * cscrr;
 yi[n - 1] = s2i * cscrr;
 if (n === 1) break mainExecutionLoop;
 nl = n - 1;
 fnui = nl;
 k = nl;
 // do 40 i=1,nl
 for (i = 1; i <= nl; i++) {
 str = s2r;
 sti = s2i;
 s2r = (fnu + fnui) * (rzr * str - rzi * sti) + s1r;
 s2i = (fnu + fnui) * (rzr * sti + rzi * str) + s1i;
 s1r = str;
 s1i = sti;
 str = s2r * cscrr;
 sti = s2i * cscrr;
 yr[k - 1] = str;
 yi[k - 1] = sti;
 fnui = fnui - 1.0;
 k = k - 1;
 if (iflag >= 3) continue;
 c1r = Math.abs(str);
 c1i = Math.abs(sti);
 c1m = Math.max(c1r, c1i);
 if (c1m <= ascle) continue;
 iflag = iflag + 1;
 ascle = bry[iflag - 1];
 s1r = s1r * cscrr;
 s1i = s1i * cscrr;
 s2r = str;
 s2i = sti;
 csclr = csclr * tol;
 cscrr = 1.0 / csclr;
 s1r = s1r * csclr;
 s1i = s1i * csclr;
 s2r = s2r * csclr;
 s2i = s2i * csclr;
 }
 // 40 continue
 break mainExecutionLoop;
 case 50:
 nz = -1;
 if (nw === -2) nz = -2;
 break mainExecutionLoop;
 case 60:
 if (iform === 2) {
 goToLabel = 70;break;
 }
 // c-----------------------------------------------------------------------
 // c asymptotic expansion for i(fnu,z) for large fnu applied in
 // c -pi/3 <= arg(z) <= pi/3
 // c-----------------------------------------------------------------------

 var _zuni5 = (0, _zuni9.zuni1)(zr, zi, fnu, kode, n, yr, yi, fnul, tol, elim, alim);

 var _zuni6 = _slicedToArray(_zuni5, 2);

 nw = _zuni6[0];
 nlast = _zuni6[1];

 goToLabel = 80;break;
 case 70:
 var _zuni7 = (0, _zuni10.zuni2)(zr, zi, fnu, kode, n, yr, yi, fnul, tol, elim, alim);
 // c-----------------------------------------------------------------------
 // c asymptotic expansion for j(fnu,z*exp(m*hpi)) for large fnu
 // c applied in pi/3 < abs(arg(z)) <= pi/2 where m=+i or -i
 // c and hpi=pi/2
 // c-----------------------------------------------------------------------


 var _zuni8 = _slicedToArray(_zuni7, 2);

 nw = _zuni8[0];
 nlast = _zuni8[1];

 case 80:
 if (nw < 0) {
 goToLabel = 50;break;
 }
 nz = nw;
 break mainExecutionLoop;
 case 90:
 nlast = n;
 default:
 break mainExecutionLoop;
 }
 }

 return [nz, nlast];
}
},{"../../utils/fortran-utils/d1mach.js":91,"./zabs.js":11,"./zuni1.js":39,"./zuni2.js":40}],25:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.zbunk = zbunk;

var _zunk = require('./zunk1.js');

var _zunk2 = require('./zunk2.js');

// SUBROUTINE ZBUNK(ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, TOL, ELIM, ALIM)
// ***BEGIN PROLOGUE ZBUNK
// ***REFER TO ZBESK,ZBESH
//
// ZBUNK COMPUTES THE K BESSEL FUNCTION FOR FNU.GT.FNUL.
// ACCORDING TO THE UNIFORM ASYMPTOTIC EXPANSION FOR K(FNU,Z)
// IN ZUNK1 AND THE EXPANSION FOR H(2,FNU,Z) IN ZUNK2
//
// ***ROUTINES CALLED ZUNK1,ZUNK2
// ***END PROLOGUE ZBUNK
function zbunk(zr, zi, fnu, kode, mr, n, yr, yi, tol, elim, alim) {
 var ax = void 0,
 ay = void 0;
 ax = Math.abs(zr) * 1.7321;
 ay = Math.abs(zi);
 if (ay > ax) {
 return (0, _zunk2.zunk2)(zr, zi, fnu, kode, mr, n, yr, yi, tol, elim, alim);
 } else {
 return (0, _zunk.zunk1)(zr, zi, fnu, kode, mr, n, yr, yi, tol, elim, alim);
 }
}
},{"./zunk1.js":42,"./zunk2.js":43}],26:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.zdiv = zdiv;

var _zabs = require('./zabs.js');

function zdiv(ar, ai, br, bi) {
 var cr = void 0,
 ci = void 0,
 bm = void 0,
 ca = void 0,
 cb = void 0,
 cc = void 0,
 cd = void 0;
 bm = 1 / (0, _zabs.azabs)(br, bi);
 cc = br * bm;
 cd = bi * bm;
 ca = (ar * cc + ai * cd) * bm;
 cb = (ai * cc - ar * cd) * bm;
 cr = ca;
 ci = cb;
 return [cr, ci];
} // SUBROUTINE ZDIV(AR, AI, BR, BI, CR, CI)
// ***BEGIN PROLOGUE ZDIV
// ***REFER TO ZBESH,ZBESI,ZBESJ,ZBESK,ZBESY,ZAIRY,ZBIRY
//
// DOUBLE PRECISION COMPLEX DIVIDE C=A/B.
//
// ***ROUTINES CALLED AZABS
// ***END PROLOGUE ZDIV
},{"./zabs.js":11}],27:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.azexp = azexp;
// SUBROUTINE AZEXP(AR, AI, BR, BI)
// ***BEGIN PROLOGUE AZEXP
// ***REFER TO ZBESH,ZBESI,ZBESJ,ZBESK,ZBESY,ZAIRY,ZBIRY
//
// DOUBLE PRECISION COMPLEX EXPONENTIAL FUNCTION B=EXP(A)
//
// ***ROUTINES CALLED (NONE)
// ***END PROLOGUE AZEXP
function azexp(ar, ai) {
 var zm = Math.exp(ar);
 return [zm * Math.cos(ai), zm * Math.sin(ai)];
}
},{}],28:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});

var _slicedToArray = function () { function sliceIterator(arr, i) { var _arr = []; var _n = true; var _d = false; var _e = undefined; try { for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) { _arr.push(_s.value); if (i && _arr.length === i) break; } } catch (err) { _d = true; _e = err; } finally { try { if (!_n && _i["return"]) _i["return"](); } finally { if (_d) throw _e; } } return _arr; } return function (arr, i) { if (Array.isArray(arr)) { return arr; } else if (Symbol.iterator in Object(arr)) { return sliceIterator(arr, i); } else { throw new TypeError("Invalid attempt to destructure non-iterable instance"); } }; }(); /* eslint-disable no-labels */
/* eslint-disable no-fallthrough */
//* **BEGIN PROLOGUE ZKSCL
//* **REFER TO ZBESK
//
// SET K FUNCTIONS TO ZERO ON UNDERFLOW, RECURRENCE
// ON SCALED FUNCTIONS UNTIL TWO MEMBERS COME ON SCALE, {
// RETURN WITH MIN(NZ+2,N) VALUES SCALED BY 1/TOL.
//
//* **ROUTINES ED() ZUCHK,AZABS,AZLOG
//* **END PROLOGUE ZKSCL
// COMPLEX CK,CS,CY,CZERO,RZ,S1,S2,Y,ZR,ZD,CELM


exports.zkscl = zkscl;

var _zabs = require('./zabs.js');

var _zlog = require('./zlog.js');

var _zuchk = require('./zuchk.js');

function zkscl(zrr, zri, fnu, n, yr, yi, nz, rzr, rzi, ascle, tol, elim) {
 var acs = void 0,
 as = void 0,
 cki = void 0,
 ckr = void 0,
 csi = void 0,
 csr = void 0,
 cyi = void 0,
 cyr = void 0,
 fn = void 0,
 str = void 0,
 s1i = void 0,
 s1r = void 0,
 s2i = void 0,
 s2r = void 0,
 zeroi = void 0,
 zeror = void 0,
 zdr = void 0,
 zdi = void 0,
 celmr = void 0,
 elm = void 0,
 helim = void 0,
 alas = void 0,
 i = void 0,
 ic = void 0,
 kk = void 0,
 nn = void 0,
 nw = void 0;

 cyr = new Float64Array(2);
 cyi = new Float64Array(2);

 var goToLabel = 0;
 // console.log('zkscl');
 mainExecutionLoop: while (true) {
 switch (goToLabel) {
 case 0:
 //
 zeror = 0.0e0;
 // double precision acs, as, ascle, cki, ckr, csi, csr, cyi, cyr, elim, fn, fnu, rzi, rzr, str, s1i, s1r, s2i, s2r, tol, yi, yr, zeroi, zeror, zri, zrr, azabs, zdr, zdi, celmr, elm, helim, alas
 // integer i, ic, kk, n, nn, nw, nz
 // dimension yr(n), yi(n), cyr(2), cyi(2)

 zeroi = 0.0e0;
 nz = 0;
 ic = 0;
 nn = Math.min(2, n);
 for (i = 1; i <= nn; i++) {
 s1r = yr[i - 1];
 s1i = yi[i - 1];
 cyr[i - 1] = s1r;
 cyi[i - 1] = s1i;
 as = (0, _zabs.azabs)(s1r, s1i);
 acs = -zrr + Math.log(as);
 nz = nz + 1;
 yr[i - 1] = zeror;
 yi[i - 1] = zeroi;
 if (acs < -elim) continue;

 var _azlog = (0, _zlog.azlog)(s1r, s1i);

 var _azlog2 = _slicedToArray(_azlog, 2);

 csr = _azlog2[0];
 csi = _azlog2[1];

 csr = csr - zrr;
 csi = csi - zri;
 str = Math.exp(csr) / tol;
 csr = str * Math.cos(csi);
 csi = str * Math.sin(csi);
 (0, _zuchk.zuchk)(csr, csi, nw, ascle, tol);
 if (nw !== 0) continue;
 yr[i - 1] = csr;
 yi[i - 1] = csi;
 ic = i;
 nz = nz - 1;
 }
 if (n === 1) break mainExecutionLoop;
 if (ic > 1) {
 goToLabel = 20;break;
 }
 yr[0] = zeror;
 yi[0] = zeroi;
 nz = 2;
 case 20:
 if (n === 2) break mainExecutionLoop;
 if (nz === 0) break mainExecutionLoop;
 fn = fnu + 1.0e0;
 ckr = fn * rzr;
 cki = fn * rzi;
 s1r = cyr[0];
 s1i = cyi[0];
 s2r = cyr[1];
 s2i = cyi[1];
 helim = 0.5e0 * elim;
 elm = Math.exp(-elim);
 celmr = elm;
 zdr = zrr;
 zdi = zri;
 //
 // find two consecutive y values on scale. scale recurrence if
 // s2 gets larger than Math.exp(elim/2)
 //
 for (i = 3; i <= n; i++) {
 kk = i;
 csr = s2r;
 csi = s2i;
 s2r = ckr * csr - cki * csi + s1r;
 s2i = cki * csr + ckr * csi + s1i;
 s1r = csr;
 s1i = csi;
 ckr = ckr + rzr;
 cki = cki + rzi;
 as = (0, _zabs.azabs)(s2r, s2i);
 alas = Math.log(as);
 acs = -zdr + alas;
 nz = nz + 1;
 yr[i - 1] = zeror;
 yi[i - 1] = zeroi;
 if (acs < -elim) {
 // goToLabel = 25; break;
 } else {
 var _azlog3 = (0, _zlog.azlog)(s2r, s2i);

 var _azlog4 = _slicedToArray(_azlog3, 2);

 csr = _azlog4[0];
 csi = _azlog4[1];

 csr = csr - zdr;
 csi = csi - zdi;
 str = Math.exp(csr) / tol;
 csr = str * Math.cos(csi);
 csi = str * Math.sin(csi);
 (0, _zuchk.zuchk)(csr, csi, nw, ascle, tol);
 if (nw !== 0) {
 // goToLabel = 25; break;
 } else {
 yr[i - 1] = csr;
 yi[i - 1] = csi;
 nz = nz - 1;
 if (ic === kk - 1) {
 goToLabel = 40;break;
 }
 ic = kk;
 continue;
 }
 }
 // case 25:
 if (alas < helim) continue;
 zdr = zdr - elim;
 s1r = s1r * celmr;
 s1i = s1i * celmr;
 s2r = s2r * celmr;
 s2i = s2i * celmr;
 }
 if (goToLabel > 30) break;
 nz = n;
 if (ic === n) nz = n - 1;
 goToLabel = 45;break;
 case 40:
 nz = kk - 2;
 case 45:
 for (i = 1; i <= nz; i++) {
 yr[i - 1] = zeror;
 yi[i - 1] = zeroi;
 }
 default:
 break mainExecutionLoop;
 }
 }
 return nz;
}
},{"./zabs.js":11,"./zlog.js":29,"./zuchk.js":37}],29:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.azlog = azlog;

var _zabs = require('./zabs.js');

function azlog(ar, ai, br, bi, ierr) {
 var zm = void 0,
 dtheta = void 0,
 dpi = void 0,
 dhpi = void 0;
 var goToLabel = 0;
 mainExecutionLoop: while (true) {
 switch (goToLabel) {
 case 0:
 dpi = 3.141592653589793238462643383;
 // double precision ar, ai, br, bi, zm, dtheta, dpi, dhpi
 // double precision azabs

 dhpi = 1.570796326794896619231321696;

 ierr = 0;
 if (ar === 0.0e+0) {
 goToLabel = 10;break;
 }
 if (ai === 0.0e+0) {
 goToLabel = 20;break;
 }
 dtheta = Math.atan(ai / ar);
 if (dtheta <= 0.0e+0) {
 goToLabel = 40;break;
 }
 if (ar < 0.0e+0) dtheta = dtheta - dpi;
 goToLabel = 50;break;
 case 10:
 if (ai === 0.0e+0) {
 goToLabel = 60;break;
 }
 bi = dhpi;
 br = Math.log(Math.abs(ai));
 if (ai < 0.0e+0) bi = -bi;
 break mainExecutionLoop;
 case 20:
 if (ar > 0.0e+0) {
 goToLabel = 30;break;
 }
 br = Math.log(Math.abs(ar));
 bi = dpi;
 break mainExecutionLoop;
 case 30:
 br = Math.log(ar);
 bi = 0.0e+0;
 break mainExecutionLoop;
 case 40:
 if (ar < 0.0e+0) dtheta = dtheta + dpi;
 case 50:
 zm = (0, _zabs.azabs)(ar, ai);
 br = Math.log(zm);
 bi = dtheta;
 break mainExecutionLoop;
 case 60:

 ierr = 1;
 default:
 break mainExecutionLoop;
 }
 }
 return [br, bi, ierr];
} /* eslint-disable no-labels */
/* eslint-disable no-fallthrough */
// DOUBLE PRECISION COMPLEX LOGARITHM B=CLOG(A)
// IERR=0,NORMAL RETURN IERR=1, Z=CMPLX(0.0,0.0)
},{"./zabs.js":11}],30:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});

var _slicedToArray = function () { function sliceIterator(arr, i) { var _arr = []; var _n = true; var _d = false; var _e = undefined; try { for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) { _arr.push(_s.value); if (i && _arr.length === i) break; } } catch (err) { _d = true; _e = err; } finally { try { if (!_n && _i["return"]) _i["return"](); } finally { if (_d) throw _e; } } return _arr; } return function (arr, i) { if (Array.isArray(arr)) { return arr; } else if (Symbol.iterator in Object(arr)) { return sliceIterator(arr, i); } else { throw new TypeError("Invalid attempt to destructure non-iterable instance"); } }; }(); /* eslint-disable no-labels */
/* eslint-disable no-fallthrough */
// SUBROUTINE ZMLRI(ZR, ZI, FNU, KODE, N, YR, YI, NZ, TOL)
// ***BEGIN PROLOGUE ZMLRI
// ***REFER TO ZBESI,ZBESK
//
// ZMLRI COMPUTES THE I BESSEL FUNCTION FOR RE(Z).GE.0.0 BY THE
// MILLER ALGORITHM NORMALIZED BY A NEUMANN SERIES.
//
// ***ROUTINES CALLED DGAMLN,D1MACH,AZABS,AZEXP,AZLOG,ZMLT
// ***END PROLOGUE ZMLRI


exports.zmlri = zmlri;

var _dgamln = require('./dgamln.js');

var _d1mach = require('../../utils/fortran-utils/d1mach.js');

var _zabs = require('./zabs.js');

var _zexp = require('./zexp.js');

var _zlog = require('./zlog.js');

var _zmlt3 = require('./zmlt.js');

function zmlri(zr, zi, fnu, kode, n, yr, yi, tol) {
 var ack = void 0,
 ak = void 0,
 ap = void 0,
 at = void 0,
 az = void 0,
 bk = void 0,
 cki = void 0,
 ckr = void 0,
 cnormi = void 0,
 cnormr = void 0,
 conei = void 0,
 coner = void 0,
 fkap = void 0,
 fkk = void 0,
 flam = void 0,
 fnf = void 0,
 pti = void 0,
 ptr = void 0,
 p1i = void 0,
 p1r = void 0,
 p2i = void 0,
 p2r = void 0,
 raz = void 0,
 rho = void 0,
 rho2 = void 0,
 rzi = void 0,
 rzr = void 0,
 scle = void 0,
 sti = void 0,
 str = void 0,
 sumi = void 0,
 sumr = void 0,
 tfnf = void 0,
 tst = void 0,
 zeroi = void 0,
 zeror = void 0,
 i = void 0,
 iaz = void 0,
 ifnu = void 0,
 inu = void 0,
 itime = void 0,
 k = void 0,
 kk = void 0,
 km = void 0,
 m = void 0,
 nz = void 0;
 zeror = 0.0;
 zeroi = 0.0;
 coner = 1.0;
 conei = 0.0;


 var goToLabel = 0;
 mainExecutionLoop: while (true) {
 mainSwitch: switch (goToLabel) {
 case 0:
 scle = (0, _d1mach.d1mach)(1) / tol;
 nz = 0;
 az = (0, _zabs.azabs)(zr, zi);
 iaz = Math.trunc(az);
 ifnu = Math.trunc(fnu);
 inu = ifnu + n - 1;
 at = iaz + 1.0;
 raz = 1.0 / az;
 str = zr * raz;
 sti = -zi * raz;
 ckr = str * at * raz;
 cki = sti * at * raz;
 rzr = (str + str) * raz;
 rzi = (sti + sti) * raz;
 p1r = zeror;
 p1i = zeroi;
 p2r = coner;
 p2i = conei;
 ack = (at + 1.0) * raz;
 rho = ack + Math.sqrt(ack * ack - 1.0);
 rho2 = rho * rho;
 tst = (rho2 + rho2) / ((rho2 - 1.0) * (rho - 1.0));
 tst = tst / tol;
 // c-----------------------------------------------------------------------
 // c compute relative truncation error index for series
 // c-----------------------------------------------------------------------
 ak = at;
 // do 10 i=1,80
 for (i = 1; i <= 80; i++) {
 ptr = p2r;
 pti = p2i;
 p2r = p1r - (ckr * ptr - cki * pti);
 p2i = p1i - (cki * ptr + ckr * pti);
 p1r = ptr;
 p1i = pti;
 ckr = ckr + rzr;
 cki = cki + rzi;
 ap = (0, _zabs.azabs)(p2r, p2i);
 if (ap > tst * ak * ak) {
 goToLabel = 20;break mainSwitch;
 }
 ak = ak + 1.0;
 }
 goToLabel = 110;break;
 case 20:
 i = i + 1;
 k = 0;
 if (inu < iaz) {
 goToLabel = 40;break;
 }
 // c-----------------------------------------------------------------------
 // c compute relative truncation error for ratios
 // c-----------------------------------------------------------------------
 p1r = zeror;
 p1i = zeroi;
 p2r = coner;
 p2i = conei;
 at = inu + 1.0;
 str = zr * raz;
 sti = -zi * raz;
 ckr = str * at * raz;
 cki = sti * at * raz;
 ack = at * raz;
 tst = Math.sqrt(ack / tol);
 itime = 1;
 // do 30 k=1,80
 for (k = 1; k <= 80; k++) {
 ptr = p2r;
 pti = p2i;
 p2r = p1r - (ckr * ptr - cki * pti);
 p2i = p1i - (ckr * pti + cki * ptr);
 p1r = ptr;
 p1i = pti;
 ckr = ckr + rzr;
 cki = cki + rzi;
 ap = (0, _zabs.azabs)(p2r, p2i);
 if (ap < tst) continue;
 if (itime === 2) {
 goToLabel = 40;break;
 }
 ack = (0, _zabs.azabs)(ckr, cki);
 flam = ack + Math.sqrt(ack * ack - 1.0);
 fkap = ap / (0, _zabs.azabs)(p1r, p1i);
 rho = Math.min(flam, fkap);
 tst = tst * Math.sqrt(rho / (rho * rho - 1.0));
 itime = 2;
 }
 if (goToLabel < 40) {
 goToLabel = 110;break;
 }
 case 40:
 // c-----------------------------------------------------------------------
 // c backward recurrence and sum normalizing relation
 // c-----------------------------------------------------------------------
 k = k + 1;
 kk = Math.max(i + iaz, k + inu);
 fkk = kk;
 p1r = zeror;
 p1i = zeroi;
 // c-----------------------------------------------------------------------
 // c scale p2 and sum by scle
 // c-----------------------------------------------------------------------
 p2r = scle;
 p2i = zeroi;
 fnf = fnu - ifnu;
 tfnf = fnf + fnf;
 bk = (0, _dgamln.dgamln)(fkk + tfnf + 1.0) - (0, _dgamln.dgamln)(fkk + 1.0) - (0, _dgamln.dgamln)(tfnf + 1.0);
 bk = Math.exp(bk);
 sumr = zeror;
 sumi = zeroi;
 km = kk - inu;
 // do 50 i=1,km
 for (i = 1; i <= km; i++) {
 ptr = p2r;
 pti = p2i;
 p2r = p1r + (fkk + fnf) * (rzr * ptr - rzi * pti);
 p2i = p1i + (fkk + fnf) * (rzi * ptr + rzr * pti);
 p1r = ptr;
 p1i = pti;
 ak = 1.0 - tfnf / (fkk + tfnf);
 ack = bk * ak;
 sumr = sumr + (ack + bk) * p1r;
 sumi = sumi + (ack + bk) * p1i;
 bk = ack;
 fkk = fkk - 1.0;
 }
 yr[n - 1] = p2r;
 yi[n - 1] = p2i;
 if (n === 1) {
 goToLabel = 70;break;
 }
 // do 60 i=2,n
 for (i = 2; i <= n; i++) {
 ptr = p2r;
 pti = p2i;
 p2r = p1r + (fkk + fnf) * (rzr * ptr - rzi * pti);
 p2i = p1i + (fkk + fnf) * (rzi * ptr + rzr * pti);
 p1r = ptr;
 p1i = pti;
 ak = 1.0 - tfnf / (fkk + tfnf);
 ack = bk * ak;
 sumr = sumr + (ack + bk) * p1r;
 sumi = sumi + (ack + bk) * p1i;
 bk = ack;
 fkk = fkk - 1.0;
 m = n - i + 1;
 yr[m - 1] = p2r;
 yi[m - 1] = p2i;
 }
 case 70:
 if (ifnu <= 0) {
 goToLabel = 90;break;
 }
 // do 80 i=1,ifnu
 for (i = 1; i <= fnu; i++) {
 ptr = p2r;
 pti = p2i;
 p2r = p1r + (fkk + fnf) * (rzr * ptr - rzi * pti);
 p2i = p1i + (fkk + fnf) * (rzr * pti + rzi * ptr);
 p1r = ptr;
 p1i = pti;
 ak = 1.0 - tfnf / (fkk + tfnf);
 ack = bk * ak;
 sumr = sumr + (ack + bk) * p1r;
 sumi = sumi + (ack + bk) * p1i;
 bk = ack;
 fkk = fkk - 1.0;
 }
 case 90:
 ptr = zr;
 pti = zi;
 if (kode === 2) ptr = zeror;

 var _azlog = (0, _zlog.azlog)(rzr, rzi);

 var _azlog2 = _slicedToArray(_azlog, 2);

 str = _azlog2[0];
 sti = _azlog2[1];

 p1r = -fnf * str + ptr;
 p1i = -fnf * sti + pti;
 ap = (0, _dgamln.dgamln)(1.0 + fnf);
 ptr = p1r - ap;
 pti = p1i;
 // c-----------------------------------------------------------------------
 // c the division cexp(pt)/(sum+p2) is altered to avoid overflow
 // c in the denominator by squaring large quantities
 // c-----------------------------------------------------------------------
 p2r = p2r + sumr;
 p2i = p2i + sumi;
 ap = (0, _zabs.azabs)(p2r, p2i);
 p1r = 1.0 / ap;

 var _azexp = (0, _zexp.azexp)(ptr, pti);

 var _azexp2 = _slicedToArray(_azexp, 2);

 str = _azexp2[0];
 sti = _azexp2[1];

 ckr = str * p1r;
 cki = sti * p1r;
 ptr = p2r * p1r;
 pti = -p2i * p1r;

 // do 100 i=1,n
 var _zmlt = (0, _zmlt3.zmlt)(ckr, cki, ptr, pti);

 var _zmlt2 = _slicedToArray(_zmlt, 2);

 cnormr = _zmlt2[0];
 cnormi = _zmlt2[1];
 for (i = 1; i <= n; i++) {
 str = yr[i - 1] * cnormr - yi[i - 1] * cnormi;
 yi[i - 1] = yr[i - 1] * cnormi + yi[i - 1] * cnormr;
 yr[i - 1] = str;
 }
 break mainExecutionLoop;
 case 110:
 nz = -2;
 default:
 break mainExecutionLoop;
 }
 }

 return nz;
}
},{"../../utils/fortran-utils/d1mach.js":91,"./dgamln.js":10,"./zabs.js":11,"./zexp.js":27,"./zlog.js":29,"./zmlt.js":31}],31:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.zmlt = zmlt;
// SUBROUTINE ZMLT(AR, AI, BR, BI, CR, CI)
// ***BEGIN PROLOGUE ZMLT
// ***REFER TO ZBESH,ZBESI,ZBESJ,ZBESK,ZBESY,ZAIRY,ZBIRY
//
// DOUBLE PRECISION COMPLEX MULTIPLY, C=A*B.
//
// ***ROUTINES CALLED (NONE)
// ***END PROLOGUE ZMLT
function zmlt(ar, ai, br, bi) {
 return [ar * br - ai * bi, ar * bi + ai * br];
}
},{}],32:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});

var _slicedToArray = function () { function sliceIterator(arr, i) { var _arr = []; var _n = true; var _d = false; var _e = undefined; try { for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) { _arr.push(_s.value); if (i && _arr.length === i) break; } } catch (err) { _d = true; _e = err; } finally { try { if (!_n && _i["return"]) _i["return"](); } finally { if (_d) throw _e; } } return _arr; } return function (arr, i) { if (Array.isArray(arr)) { return arr; } else if (Symbol.iterator in Object(arr)) { return sliceIterator(arr, i); } else { throw new TypeError("Invalid attempt to destructure non-iterable instance"); } }; }(); /* eslint-disable no-labels */
/* eslint-disable no-fallthrough */
// SUBROUTINE ZRATI(ZR, ZI, FNU, N, CYR, CYI, TOL)
// ***BEGIN PROLOGUE ZRATI
// ***REFER TO ZBESI,ZBESK,ZBESH
//
// ZRATI COMPUTES RATIOS OF I BESSEL FUNCTIONS BY BACKWARD
// RECURRENCE. THE STARTING INDEX IS DETERMINED BY FORWARD
// RECURRENCE AS DESCRIBED IN J. RES. OF NAT. BUR. OF STANDARDS-B,
// MATHEMATICAL SCIENCES, VOL 77B, P111-114, SEPTEMBER, 1973,
// BESSEL FUNCTIONS I AND J OF COMPLEX ARGUMENT AND INTEGER ORDER,
// BY D. J. SOOKNE.
//
// ***ROUTINES CALLED AZABS,ZDIV
// ***END PROLOGUE ZRATI


exports.zrati = zrati;

var _zabs = require('./zabs.js');

var _zdiv3 = require('./zdiv.js');

function zrati(zr, zi, fnu, n, cyr, cyi, tol) {
 var ak = void 0,
 amagz = void 0,
 ap1 = void 0,
 ap2 = void 0,
 arg = void 0,
 az = void 0,
 cdfnui = void 0,
 cdfnur = void 0,
 conei = void 0,
 coner = void 0,
 czeroi = void 0,
 czeror = void 0,
 dfnu = void 0,
 fdnu = void 0,
 flam = void 0,
 fnup = void 0,
 pti = void 0,
 ptr = void 0,
 p1i = void 0,
 p1r = void 0,
 p2i = void 0,
 p2r = void 0,
 rak = void 0,
 rap1 = void 0,
 rho = void 0,
 rt2 = void 0,
 rzi = void 0,
 rzr = void 0,
 test = void 0,
 test1 = void 0,
 tti = void 0,
 ttr = void 0,
 t1i = void 0,
 t1r = void 0,
 i = void 0,
 id = void 0,
 idnu = void 0,
 inu = void 0,
 itime = void 0,
 k = void 0,
 kk = void 0,
 magz = void 0;

 czeror = 0.0;
 czeroi = 0.0;
 coner = 1.0;
 conei = 0.0;
 rt2 = 1.41421356237309505;


 var goToLabel = 0;
 mainExecutionLoop: while (true) {
 switch (goToLabel) {
 case 0:
 az = (0, _zabs.azabs)(zr, zi);
 inu = Math.trunc(fnu);
 idnu = inu + n - 1;
 magz = Math.trunc(az);
 amagz = magz + 1;
 fdnu = idnu;
 fnup = Math.max(amagz, fdnu);
 id = idnu - magz - 1;
 itime = 1;
 k = 1;
 ptr = 1.0 / az;
 rzr = ptr * (zr + zr) * ptr;
 rzi = -ptr * (zi + zi) * ptr;
 t1r = rzr * fnup;
 t1i = rzi * fnup;
 p2r = -t1r;
 p2i = -t1i;
 p1r = coner;
 p1i = conei;
 t1r = t1r + rzr;
 t1i = t1i + rzi;
 if (id > 0) id = 0;
 ap2 = (0, _zabs.azabs)(p2r, p2i);
 ap1 = (0, _zabs.azabs)(p1r, p1i);
 // c-----------------------------------------------------------------------
 // c the overflow test on k(fnu+i-1,z) before the call to cbknu
 // c guarantees that p2 is on scale. scale test1 and all subsequent
 // c p2 values by ap1 to ensure that an overflow does not occur
 // c prematurely.
 // c-----------------------------------------------------------------------
 arg = (ap2 + ap2) / (ap1 * tol);
 test1 = Math.sqrt(arg);
 test = test1;
 rap1 = 1.0 / ap1;
 p1r = p1r * rap1;
 p1i = p1i * rap1;
 p2r = p2r * rap1;
 p2i = p2i * rap1;
 ap2 = ap2 * rap1;
 case 10:
 k = k + 1;
 ap1 = ap2;
 ptr = p2r;
 pti = p2i;
 p2r = p1r - (t1r * ptr - t1i * pti);
 p2i = p1i - (t1r * pti + t1i * ptr);
 p1r = ptr;
 p1i = pti;
 t1r = t1r + rzr;
 t1i = t1i + rzi;
 ap2 = (0, _zabs.azabs)(p2r, p2i);
 if (ap1 <= test) {
 goToLabel = 10;break;
 }
 if (itime === 2) {
 goToLabel = 20;break;
 }
 ak = (0, _zabs.azabs)(t1r, t1i) * 0.5;
 flam = ak + Math.sqrt(ak * ak - 1.0);
 rho = Math.min(ap2 / ap1, flam);
 test = test1 * Math.sqrt(rho / (rho * rho - 1.0));
 itime = 2;
 goToLabel = 10;break;
 case 20:
 kk = k + 1 - id;
 ak = kk;
 t1r = ak;
 t1i = czeroi;
 dfnu = fnu + (n - 1);
 p1r = 1.0 / ap2;
 p1i = czeroi;
 p2r = czeror;
 p2i = czeroi;
 // do 30 i=1,kk
 for (i = 1; i <= kk; i++) {
 ptr = p1r;
 pti = p1i;
 rap1 = dfnu + t1r;
 ttr = rzr * rap1;
 tti = rzi * rap1;
 p1r = ptr * ttr - pti * tti + p2r;
 p1i = ptr * tti + pti * ttr + p2i;
 p2r = ptr;
 p2i = pti;
 t1r = t1r - coner;
 }
 // 30 continue
 if (p1r !== czeror || p1i !== czeroi) {
 goToLabel = 40;break;
 }
 p1r = tol;
 p1i = tol;
 case 40:
 var _zdiv = (0, _zdiv3.zdiv)(p2r, p2i, p1r, p1i);

 var _zdiv2 = _slicedToArray(_zdiv, 2);

 cyr[n - 1] = _zdiv2[0];
 cyi[n - 1] = _zdiv2[1];

 if (n === 1) break mainExecutionLoop;
 k = n - 1;
 ak = k;
 t1r = ak;
 t1i = czeroi;
 cdfnur = fnu * rzr;
 cdfnui = fnu * rzi;
 // do 60 i=2,n
 for (i = 2; i <= n; i++) {
 ptr = cdfnur + (t1r * rzr - t1i * rzi) + cyr[k];
 pti = cdfnui + (t1r * rzi + t1i * rzr) + cyi[k];
 ak = (0, _zabs.azabs)(ptr, pti);
 if (ak !== czeror) {
 // go to 50
 } else {
 ptr = tol;
 pti = tol;
 ak = tol * rt2;
 }
 // 50 continue
 rak = coner / ak;
 cyr[k - 1] = rak * ptr * rak;
 cyi[k - 1] = -rak * pti * rak;
 t1r = t1r - coner;
 k = k - 1;
 }
 // 60 continue
 default:
 break mainExecutionLoop;
 }
 }
}
},{"./zabs.js":11,"./zdiv.js":26}],33:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});

var _slicedToArray = function () { function sliceIterator(arr, i) { var _arr = []; var _n = true; var _d = false; var _e = undefined; try { for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) { _arr.push(_s.value); if (i && _arr.length === i) break; } } catch (err) { _d = true; _e = err; } finally { try { if (!_n && _i["return"]) _i["return"](); } finally { if (_d) throw _e; } } return _arr; } return function (arr, i) { if (Array.isArray(arr)) { return arr; } else if (Symbol.iterator in Object(arr)) { return sliceIterator(arr, i); } else { throw new TypeError("Invalid attempt to destructure non-iterable instance"); } }; }(); // SUBROUTINE ZS1S2(ZRR, ZRI, S1R, S1I, S2R, S2I, NZ, ASCLE, ALIM, IUF)
// ***BEGIN PROLOGUE ZS1S2
// ***REFER TO ZBESK,ZAIRY
//
// ZS1S2 TESTS FOR A POSSIBLE UNDERFLOW RESULTING FROM THE
// ADDITION OF THE I AND K FUNCTIONS IN THE ANALYTIC CON-
// TINUATION FORMULA WHERE S1=K FUNCTION AND S2=I FUNCTION.
// ON KODE=1 THE I AND K FUNCTIONS ARE DIFFERENT ORDERS OF
// MAGNITUDE, BUT FOR KODE=2 THEY CAN BE OF THE SAME ORDER
// OF MAGNITUDE AND THE MAXIMUM MUST BE AT LEAST ONE
// PRECISION ABOVE THE UNDERFLOW LIMIT.
//
// ***ROUTINES CALLED AZABS,AZEXP,AZLOG
// ***END PROLOGUE ZS1S2


exports.zs1s2 = zs1s2;

var _zabs = require('./zabs.js');

var _zexp = require('./zexp.js');

var _zlog = require('./zlog.js');

function zs1s2(zrr, zri, s1r, s1i, s2r, s2i, ascle, alim, iuf) {
 var aa = void 0,
 aln = void 0,
 as1 = void 0,
 as2 = void 0,
 c1i = void 0,
 c1r = void 0,
 s1di = void 0,
 s1dr = void 0,
 zeroi = void 0,
 zeror = void 0,
 nz = void 0;
 zeror = 0;
 zeroi = 0;
 nz = 0;
 as1 = (0, _zabs.azabs)(s1r, s1i);
 as2 = (0, _zabs.azabs)(s2r, s2i);
 if (s1r === 0 && s1i === 0 || as1 === 0) {
 // go to 10
 } else {
 aln = -zrr - zrr + Math.log(as1);
 s1dr = s1r;
 s1di = s1i;
 s1r = zeror;
 s1i = zeroi;
 as1 = zeror;
 if (aln < -alim) {
 // go to 10
 } else {
 var _azlog = (0, _zlog.azlog)(s1dr, s1di);

 var _azlog2 = _slicedToArray(_azlog, 2);

 c1r = _azlog2[0];
 c1i = _azlog2[1];

 c1r = c1r - zrr - zrr;
 c1i = c1i - zri - zri;

 var _azexp = (0, _zexp.azexp)(c1r, c1i);

 var _azexp2 = _slicedToArray(_azexp, 2);

 s1r = _azexp2[0];
 s1i = _azexp2[1];

 as1 = (0, _zabs.azabs)(s1r, s1i);
 iuf = iuf + 1;
 }
 }
 // 10 continue
 aa = Math.max(as1, as2);
 if (aa > ascle) {
 return [s1r, s1i, s2r, s2i, nz, iuf];
 }
 s1r = zeror;
 s1i = zeroi;
 s2r = zeror;
 s2i = zeroi;
 nz = 1;
 iuf = 0;
 return [s1r, s1i, s2r, s2i, nz, iuf];
}
},{"./zabs.js":11,"./zexp.js":27,"./zlog.js":29}],34:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});

var _slicedToArray = function () { function sliceIterator(arr, i) { var _arr = []; var _n = true; var _d = false; var _e = undefined; try { for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) { _arr.push(_s.value); if (i && _arr.length === i) break; } } catch (err) { _d = true; _e = err; } finally { try { if (!_n && _i["return"]) _i["return"](); } finally { if (_d) throw _e; } } return _arr; } return function (arr, i) { if (Array.isArray(arr)) { return arr; } else if (Symbol.iterator in Object(arr)) { return sliceIterator(arr, i); } else { throw new TypeError("Invalid attempt to destructure non-iterable instance"); } }; }(); /* eslint-disable no-labels */
/* eslint-disable no-fallthrough */
// SUBROUTINE ZSERI(ZR, ZI, FNU, KODE, N, YR, YI, NZ, TOL, ELIM, ALIM)
// ***BEGIN PROLOGUE ZSERI
// ***REFER TO ZBESI,ZBESK
//
// ZSERI COMPUTES THE I BESSEL FUNCTION FOR REAL(Z).GE.0.0 BY
// MEANS OF THE POWER SERIES FOR LARGE CABS(Z) IN THE
// REGION CABS(Z).LE.2*SQRT(FNU+1). NZ=0 IS A NORMAL RETURN.
// NZ.GT.0 MEANS THAT THE LAST NZ COMPONENTS WERE SET TO ZERO
// DUE TO UNDERFLOW. NZ.LT.0 MEANS UNDERFLOW OCCURRED, BUT THE
// CONDITION CABS(Z).LE.2*SQRT(FNU+1) WAS VIOLATED AND THE
// COMPUTATION MUST BE COMPLETED IN ANOTHER ROUTINE WITH N=N-ABS(NZ).
//
// ***ROUTINES CALLED DGAMLN,D1MACH,ZUCHK,AZABS,ZDIV,AZLOG,ZMLT
// ***END PROLOGUE ZSERI


exports.zseri = zseri;

var _dgamln = require('./dgamln.js');

var _d1mach = require('../../utils/fortran-utils/d1mach.js');

var _zuchk = require('./zuchk.js');

var _zabs = require('./zabs.js');

var _zdiv3 = require('./zdiv.js');

var _zlog = require('./zlog.js');

var _zmlt3 = require('./zmlt.js');

function zseri(zr, zi, fnu, kode, n, yr, yi, tol, elim, alim) {
 var aa = void 0,
 acz = void 0,
 ak = void 0,
 ak1i = void 0,
 ak1r = void 0,
 arm = void 0,
 ascle = void 0,
 atol = void 0,
 az = void 0,
 cki = void 0,
 ckr = void 0,
 coefi = void 0,
 coefr = void 0,
 conei = void 0,
 coner = void 0,
 crscr = void 0,
 czi = void 0,
 czr = void 0,
 dfnu = void 0,
 fnup = void 0,
 hzi = void 0,
 hzr = void 0,
 raz = void 0,
 rs = void 0,
 rtr1 = void 0,
 rzi = void 0,
 rzr = void 0,
 s = void 0,
 ss = void 0,
 sti = void 0,
 str = void 0,
 s1i = void 0,
 s1r = void 0,
 s2i = void 0,
 s2r = void 0,
 wi = void 0,
 wr = void 0,
 zeroi = void 0,
 zeror = void 0,
 i = void 0,
 ib = void 0,
 iflag = void 0,
 il = void 0,
 k = void 0,
 l = void 0,
 m = void 0,
 nn = void 0,
 nz = void 0,
 nw = void 0;

 wr = new Array(2);
 wi = new Array(2);

 zeror = 0.0;
 zeroi = 0.0;
 coner = 1.0;
 conei = 0.0;


 var goToLabel = 0;
 mainExecutionLoop: while (true) {
 switch (goToLabel) {
 case 0:

 nz = 0;
 az = (0, _zabs.azabs)(zr, zi);
 if (az === 0.0) {
 goToLabel = 160;break;
 }
 arm = 1.0e+3 * (0, _d1mach.d1mach)(1);
 rtr1 = Math.sqrt(arm);
 crscr = 1.0;
 iflag = 0;
 if (az < arm) {
 goToLabel = 150;break;
 }
 hzr = 0.5 * zr;
 hzi = 0.5 * zi;
 czr = zeror;
 czi = zeroi;
 if (az <= rtr1) {
 goToLabel = 10;break;
 }

 var _zmlt = (0, _zmlt3.zmlt)(hzr, hzi, hzr, hzi);

 var _zmlt2 = _slicedToArray(_zmlt, 2);

 czr = _zmlt2[0];
 czi = _zmlt2[1];

 case 10:
 acz = (0, _zabs.azabs)(czr, czi);
 nn = n;

 var _azlog = (0, _zlog.azlog)(hzr, hzi);

 var _azlog2 = _slicedToArray(_azlog, 2);

 ckr = _azlog2[0];
 cki = _azlog2[1];

 case 20:
 dfnu = fnu + (nn - 1);
 fnup = dfnu + 1.0;
 // c-----------------------------------------------------------------------
 // c underflow test
 // c-----------------------------------------------------------------------
 ak1r = ckr * dfnu;
 ak1i = cki * dfnu;
 ak = (0, _dgamln.dgamln)(fnup);
 ak1r = ak1r - ak;
 if (kode === 2) ak1r = ak1r - zr;
 if (ak1r > -elim) {
 goToLabel = 40;break;
 }
 case 30:
 nz = nz + 1;
 yr[nn - 1] = zeror;
 yi[nn - 1] = zeroi;
 if (acz > dfnu) {
 goToLabel = 190;break;
 }
 nn = nn - 1;
 if (nn === 0) break mainExecutionLoop;
 goToLabel = 20;break;
 case 40:
 if (ak1r > -alim) {
 goToLabel = 50;break;
 }
 iflag = 1;
 ss = 1.0 / tol;
 crscr = tol;
 ascle = arm * ss;
 case 50:
 aa = Math.exp(ak1r);
 if (iflag === 1) aa = aa * ss;
 coefr = aa * Math.cos(ak1i);
 coefi = aa * Math.sin(ak1i);
 atol = tol * acz / fnup;
 il = Math.min(2, nn);
 // do 90 i=1,il
 for (i = 1; i <= il; i++) {
 dfnu = fnu + (nn - i);
 fnup = dfnu + 1.0;
 s1r = coner;
 s1i = conei;
 if (acz < tol * fnup) {
 // go to 70
 } else {
 ak1r = coner;
 ak1i = conei;
 ak = fnup + 2.0;
 s = fnup;
 aa = 2.0;
 // 60 continue
 while (aa > atol) {
 rs = 1.0 / s;
 str = ak1r * czr - ak1i * czi;
 sti = ak1r * czi + ak1i * czr;
 ak1r = str * rs;
 ak1i = sti * rs;
 s1r = s1r + ak1r;
 s1i = s1i + ak1i;
 s = s + ak;
 ak = ak + 2.0;
 aa = aa * acz * rs;
 }
 // if (aa > atol) go to 60
 }
 // 70 continue
 s2r = s1r * coefr - s1i * coefi;
 s2i = s1r * coefi + s1i * coefr;
 wr[i - 1] = s2r;
 wi[i - 1] = s2i;
 if (iflag === 0) {
 // go to 80
 } else {
 nw = (0, _zuchk.zuchk)(s2r, s2i, ascle, tol);
 if (nw !== 0) {
 goToLabel = 30;break;
 }
 }
 // 80 continue
 m = nn - i + 1;
 yr[m - 1] = s2r * crscr;
 yi[m - 1] = s2i * crscr;
 if (i === il) continue;

 var _zdiv = (0, _zdiv3.zdiv)(coefr, coefi, hzr, hzi);

 var _zdiv2 = _slicedToArray(_zdiv, 2);

 str = _zdiv2[0];
 sti = _zdiv2[1];

 coefr = str * dfnu;
 coefi = sti * dfnu;
 }
 if (goToLabel === 30) break;
 case 90:
 if (nn <= 2) break mainExecutionLoop;
 k = nn - 2;
 ak = k;
 raz = 1.0 / az;
 str = zr * raz;
 sti = -zi * raz;
 rzr = (str + str) * raz;
 rzi = (sti + sti) * raz;
 if (iflag === 1) {
 goToLabel = 120;break;
 }
 ib = 3;
 case 100:
 // do 110 i=ib,nn
 for (i = ib; i <= nn; i++) {
 yr[k - 1] = (ak + fnu) * (rzr * yr[k] - rzi * yi[k]) + yr[k + 1];
 yi[k - 1] = (ak + fnu) * (rzr * yi[k] + rzi * yr[k]) + yi[k + 1];
 ak = ak - 1.0;
 k = k - 1;
 }
 break mainExecutionLoop;
 // c-----------------------------------------------------------------------
 // c recur backward with scaled values
 // c-----------------------------------------------------------------------
 case 120:
 // c-----------------------------------------------------------------------
 // c exp(-alim)=exp(-elim)/tol=approx. one precision above the
 // c underflow limit = ascle = d1mach(1)*ss*1.0e+3
 // c-----------------------------------------------------------------------
 s1r = wr[0];
 s1i = wi[0];
 s2r = wr[1];
 s2i = wi[1];
 // do 130 l=3,nn
 for (l = 3; l <= nn; l++) {
 ckr = s2r;
 cki = s2i;
 s2r = s1r + (ak + fnu) * (rzr * ckr - rzi * cki);
 s2i = s1i + (ak + fnu) * (rzr * cki + rzi * ckr);
 s1r = ckr;
 s1i = cki;
 ckr = s2r * crscr;
 cki = s2i * crscr;
 yr[k - 1] = ckr;
 yi[k - 1] = cki;
 ak = ak - 1.0;
 k = k - 1;
 if ((0, _zabs.azabs)(ckr, cki) > ascle) {
 goToLabel = 140;break;
 }
 }
 // 130 continue
 if (goToLabel !== 140) break mainExecutionLoop;
 ib = l + 1;
 if (ib > nn) break mainExecutionLoop;
 goToLabel = 100;break;
 case 150:
 nz = n;
 if (fnu === 0.0) nz = nz - 1;
 case 160:
 yr[0] = zeror;
 yi[0] = zeroi;
 if (fnu !== 0.0) {
 goToLabel = 170;break;
 }
 yr[0] = coner;
 yi[0] = conei;
 case 170:
 if (n === 1) break mainExecutionLoop;
 // do 180 i=2,n
 for (i = 2; i <= n; i++) {
 yr[i - 1] = zeror;
 yi[i - 1] = zeroi;
 }
 // 180 continue
 break mainExecutionLoop;
 // c-----------------------------------------------------------------------
 // c return with nz < 0 if cabs(z*z/4) > fnu+n-nz-1 complete
 // c the calculation in cbinu with n=n-iabs(nz)
 // c-----------------------------------------------------------------------
 case 190:
 nz = -nz;
 default:
 break mainExecutionLoop;
 }
 }

 return nz;
}
},{"../../utils/fortran-utils/d1mach.js":91,"./dgamln.js":10,"./zabs.js":11,"./zdiv.js":26,"./zlog.js":29,"./zmlt.js":31,"./zuchk.js":37}],35:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.zshch = zshch;
// SUBROUTINE ZSHCH(ZR, ZI, CSHR, CSHI, CCHR, CCHI)
// ***BEGIN PROLOGUE ZSHCH
// ***REFER TO ZBESK,ZBESH
//
// ZSHCH COMPUTES THE COMPLEX HYPERBOLIC FUNCTIONS CSH=SINH(X+I*Y)
// AND CCH=COSH(X+I*Y), WHERE I**2=-1.
//
// ***ROUTINES CALLED (NONE)
// ***END PROLOGUE ZSHCH
//
function zshch(zr, zi) {
 var sh = Math.sinh(zr);
 var ch = Math.cosh(zr);
 var sn = Math.sin(zi);
 var cn = Math.cos(zi);
 // [ cshr, cshi, cchr, cchi]
 return [sh * cn, ch * sn, ch * cn, sh * sn];
}
},{}],36:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.azsqrt = azsqrt;

var _zabs = require('./zabs.js');

function azsqrt(ar, ai) {
 var br = void 0,
 bi = void 0,
 zm = void 0,
 dtheta = void 0,
 dpi = void 0,
 drt = void 0;
 drt = 7.071067811865475244008443621e-1;
 dpi = 3.141592653589793238462643383;


 var goToLabel = 0;
 mainExecutionLoop: while (true) {
 switch (goToLabel) {
 case 0:
 zm = (0, _zabs.azabs)(ar, ai);
 zm = Math.sqrt(zm);
 if (ar === 0.0) {
 goToLabel = 10;break;
 }
 if (ai === 0.0) {
 goToLabel = 20;break;
 }
 dtheta = Math.atan(ai / ar);
 if (dtheta <= 0.0) {
 goToLabel = 40;break;
 }
 if (ar < 0.0) dtheta = dtheta - dpi;
 goToLabel = 50;break;
 case 10:
 if (ai > 0.0) {
 goToLabel = 60;break;
 }
 if (ai < 0.0) {
 goToLabel = 70;break;
 }
 br = 0.0;
 bi = 0.0;
 break mainExecutionLoop;
 case 20:
 if (ar > 0.0) {
 goToLabel = 30;break;
 }
 br = 0.0;
 bi = Math.sqrt(Math.abs(ar));
 break mainExecutionLoop;
 case 30:
 br = Math.sqrt(ar);
 bi = 0.0;
 break mainExecutionLoop;
 case 40:
 if (ar < 0.0) dtheta = dtheta + dpi;
 case 50:
 dtheta = dtheta * 0.5;
 br = zm * Math.cos(dtheta);
 bi = zm * Math.sin(dtheta);
 break mainExecutionLoop;
 case 60:
 br = zm * drt;
 bi = zm * drt;
 break mainExecutionLoop;
 case 70:
 br = zm * drt;
 bi = -zm * drt;
 default:
 break mainExecutionLoop;
 }
 }

 return [br, bi];
} /* eslint-disable no-labels */
/* eslint-disable no-fallthrough */
// SUBROUTINE AZSQRT(AR, AI, BR, BI)
// ***BEGIN PROLOGUE AZSQRT
// ***REFER TO ZBESH,ZBESI,ZBESJ,ZBESK,ZBESY,ZAIRY,ZBIRY
//
// DOUBLE PRECISION COMPLEX SQUARE ROOT, B=CSQRT(A)
//
// ***ROUTINES CALLED AZABS
// ***END PROLOGUE AZSQRT
},{"./zabs.js":11}],37:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.zuchk = zuchk;
/* eslint-disable no-labels */
/* eslint-disable no-fallthrough */
// SUBROUTINE ZUCHK(YR, YI, NZ, ASCLE, TOL)
// ***BEGIN PROLOGUE ZUCHK
// ***REFER TO ZSERI,ZUOIK,ZUNK1,ZUNK2,ZUNI1,ZUNI2,ZKSCL
//
// Y ENTERS AS A SCALED QUANTITY WHOSE MAGNITUDE IS GREATER THAN
// EXP(-ALIM)=ASCLE=1.0E+3*D1MACH(1)/TOL. THE TEST IS MADE TO SEE
// IF THE MAGNITUDE OF THE REAL OR IMAGINARY PART WOULD UNDERFLOW
// WHEN Y IS SCALED (BY TOL) TO ITS PROPER VALUE. Y IS ACCEPTED
// IF THE UNDERFLOW IS AT LEAST ONE PRECISION BELOW THE MAGNITUDE
// OF THE LARGEST COMPONENT; OTHERWISE THE PHASE ANGLE DOES NOT HAVE
// ABSOLUTE ACCURACY AND AN UNDERFLOW IS ASSUMED.
//
// ***ROUTINES CALLED (NONE)
// ***END PROLOGUE ZUCHK
//
// COMPLEX Y
function zuchk(yr, yi, ascle, tol) {
 var nz = void 0,
 ss = void 0,
 st = void 0,
 wr = void 0,
 wi = void 0;
 nz = 0;
 wr = Math.abs(yr);
 wi = Math.abs(yi);
 st = Math.min(wr, wi);
 if (st > ascle) {
 return nz;
 }
 ss = Math.max(wr, wi);
 st = st / tol;
 if (ss < st) {
 nz = 1;
 }
 return nz;
}
},{}],38:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});

var _slicedToArray = function () { function sliceIterator(arr, i) { var _arr = []; var _n = true; var _d = false; var _e = undefined; try { for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) { _arr.push(_s.value); if (i && _arr.length === i) break; } } catch (err) { _d = true; _e = err; } finally { try { if (!_n && _i["return"]) _i["return"](); } finally { if (_d) throw _e; } } return _arr; } return function (arr, i) { if (Array.isArray(arr)) { return arr; } else if (Symbol.iterator in Object(arr)) { return sliceIterator(arr, i); } else { throw new TypeError("Invalid attempt to destructure non-iterable instance"); } }; }(); /* eslint-disable no-labels */
/* eslint-disable no-fallthrough */
// SUBROUTINE ZUNHJ(ZR, ZI, FNU, IPMTR, TOL, PHIR, PHII, ARGR, ARGI,
// * ZETA1R, ZETA1I, ZETA2R, ZETA2I, ASUMR, ASUMI, BSUMR, BSUMI)
// ***BEGIN PROLOGUE ZUNHJ
// ***REFER TO ZBESI,ZBESK
//
// REFERENCES
// HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ AND I.A.
// STEGUN, AMS55, NATIONAL BUREAU OF STANDARDS, 1965, CHAPTER 9.
//
// ASYMPTOTICS AND SPECIAL FUNCTIONS BY F.W.J. OLVER, ACADEMIC
// PRESS, N.Y., 1974, PAGE 420
//
// ABSTRACT
// ZUNHJ COMPUTES PARAMETERS FOR BESSEL FUNCTIONS C(FNU,Z) =
// J(FNU,Z), Y(FNU,Z) OR H(I,FNU,Z) I=1,2 FOR LARGE ORDERS FNU
// BY MEANS OF THE UNIFORM ASYMPTOTIC EXPANSION
//
// C(FNU,Z)=C1*PHI*( ASUM*AIRY(ARG) + C2*BSUM*DAIRY(ARG) )
//
// FOR PROPER CHOICES OF C1, C2, AIRY AND DAIRY WHERE AIRY IS
// AN AIRY FUNCTION AND DAIRY IS ITS DERIVATIVE.
//
// (2/3)*FNU*ZETA**1.5 = ZETA1-ZETA2,
//
// ZETA1=0.5*FNU*CLOG((1+W)/(1-W)), ZETA2=FNU*W FOR SCALING
// PURPOSES IN AIRY FUNCTIONS FROM CAIRY OR CBIRY.
//
// MCONJ=SIGN OF AIMAG(Z), BUT IS AMBIGUOUS WHEN Z IS REAL AND
// MUST BE SPECIFIED. IPMTR=0 RETURNS ALL PARAMETERS. IPMTR=
// 1 COMPUTES ALL EXCEPT ASUM AND BSUM.
//
// ***ROUTINES CALLED AZABS,ZDIV,AZLOG,AZSQRT,D1MACH
// ***END PROLOGUE ZUNHJ
// COMPLEX ARG,ASUM,BSUM,CFNU,CONE,CR,CZERO,DR,P,PHI,PRZTH,PTFN,
// *RFN13,RTZTA,RZTH,SUMA,SUMB,TFN,T2,UP,W,W2,Z,ZA,ZB,ZC,ZETA,ZETA1,


exports.zunhj = zunhj;

var _zabs = require('./zabs.js');

var _zdiv9 = require('./zdiv.js');

var _zlog = require('./zlog.js');

var _zsqrt = require('./zsqrt.js');

var _d1mach = require('../../utils/fortran-utils/d1mach.js');

function zunhj(zr, zi, fnu, ipmtr, tol) {
 var alfa = void 0,
 ang = void 0,
 ap = void 0,
 ar = void 0,
 argi = void 0,
 argr = void 0,
 asumi = void 0,
 asumr = void 0,
 atol = void 0,
 aw2 = void 0,
 azth = void 0,
 beta = void 0,
 br = void 0,
 bsumi = void 0,
 bsumr = void 0,
 btol = void 0,
 c = void 0,
 conei = void 0,
 coner = void 0,
 cri = void 0,
 crr = void 0,
 dri = void 0,
 drr = void 0,
 ex1 = void 0,
 ex2 = void 0,
 fn13 = void 0,
 fn23 = void 0,
 gama = void 0,
 gpi = void 0,
 hpi = void 0,
 phii = void 0,
 phir = void 0,
 pi = void 0,
 pp = void 0,
 pr = void 0,
 przthi = void 0,
 przthr = void 0,
 ptfni = void 0,
 ptfnr = void 0,
 raw = void 0,
 raw2 = void 0,
 razth = void 0,
 rfnu = void 0,
 rfnu2 = void 0,
 rfn13 = void 0,
 rtzti = void 0,
 rtztr = void 0,
 rzthi = void 0,
 rzthr = void 0,
 sti = void 0,
 str = void 0,
 sumai = void 0,
 sumar = void 0,
 sumbi = void 0,
 sumbr = void 0,
 test = void 0,
 tfni = void 0,
 tfnr = void 0,
 thpi = void 0,
 tzai = void 0,
 tzar = void 0,
 t2i = void 0,
 t2r = void 0,
 upi = void 0,
 upr = void 0,
 wi = void 0,
 wr = void 0,
 w2i = void 0,
 w2r = void 0,
 zai = void 0,
 zar = void 0,
 zbi = void 0,
 zbr = void 0,
 zci = void 0,
 zcr = void 0,
 zeroi = void 0,
 zeror = void 0,
 zetai = void 0,
 zetar = void 0,
 zeta1i = void 0,
 zeta1r = void 0,
 zeta2i = void 0,
 zeta2r = void 0,
 zthi = void 0,
 zthr = void 0,
 ac = void 0,
 ias = void 0,
 ibs = void 0,
 is = void 0,
 j = void 0,
 jr = void 0,
 ju = void 0,
 k = void 0,
 kmax = void 0,
 kp1 = void 0,
 ks = void 0,
 l = void 0,
 lr = void 0,
 lrp1 = void 0,
 l1 = void 0,
 l2 = void 0,
 m = void 0;

 ap = new Array(30);
 pr = new Array(30);
 pi = new Array(30);
 upr = new Array(14);
 upi = new Array(14);
 crr = new Array(14);
 cri = new Array(14);
 drr = new Array(14);
 dri = new Array(14);

 ar = [1.00000000000000000e+00, 1.04166666666666667e-01, 8.35503472222222222e-02, 1.28226574556327160e-01, 2.91849026464140464e-01, 8.81627267443757652e-01, 3.32140828186276754e+00, 1.49957629868625547e+01, 7.89230130115865181e+01, 4.74451538868264323e+02, 3.20749009089066193e+03, 2.40865496408740049e+04, 1.98923119169509794e+05, 1.79190200777534383e+06];
 br = [1.00000000000000000e+00, -1.45833333333333333e-01, -9.87413194444444444e-02, -1.43312053915895062e-01, -3.17227202678413548e-01, -9.42429147957120249e-01, -3.51120304082635426e+00, -1.57272636203680451e+01, -8.22814390971859444e+01, -4.92355370523670524e+02, -3.31621856854797251e+03, -2.48276742452085896e+04, -2.04526587315129788e+05, -1.83844491706820990e+06];
 c = [1.00000000000000000e+00, -2.08333333333333333e-01, 1.25000000000000000e-01, 3.34201388888888889e-01, -4.01041666666666667e-01, 7.03125000000000000e-02, -1.02581259645061728e+00, 1.84646267361111111e+00, -8.91210937500000000e-01, 7.32421875000000000e-02, 4.66958442342624743e+00, -1.12070026162229938e+01, 8.78912353515625000e+00, -2.36408691406250000e+00, 1.12152099609375000e-01, -2.82120725582002449e+01, 8.46362176746007346e+01, -9.18182415432400174e+01, 4.25349987453884549e+01, -7.36879435947963170e+00, 2.27108001708984375e-01, 2.12570130039217123e+02, -7.65252468141181642e+02, 1.05999045252799988e+03, -6.99579627376132541e+02, 2.18190511744211590e+02, -2.64914304869515555e+01, 5.72501420974731445e-01, -1.91945766231840700e+03, 8.06172218173730938e+03, -1.35865500064341374e+04, 1.16553933368645332e+04, -5.30564697861340311e+03, 1.20090291321635246e+03, -1.08090919788394656e+02, 1.72772750258445740e+00, 2.02042913309661486e+04, -9.69805983886375135e+04, 1.92547001232531532e+05, -2.03400177280415534e+05, 1.22200464983017460e+05, -4.11926549688975513e+04, 7.10951430248936372e+03, -4.93915304773088012e+02, 6.07404200127348304e+00, -2.42919187900551333e+05, 1.31176361466297720e+06, -2.99801591853810675e+06, 3.76327129765640400e+06, -2.81356322658653411e+06, 1.26836527332162478e+06, -3.31645172484563578e+05, 4.52187689813627263e+04, -2.49983048181120962e+03, 2.43805296995560639e+01, 3.28446985307203782e+06, -1.97068191184322269e+07, 5.09526024926646422e+07, -7.41051482115326577e+07, 6.63445122747290267e+07, -3.75671766607633513e+07, 1.32887671664218183e+07, -2.78561812808645469e+06, 3.08186404612662398e+05, -1.38860897537170405e+04, 1.10017140269246738e+02, -4.93292536645099620e+07, 3.25573074185765749e+08, -9.39462359681578403e+08, 1.55359689957058006e+09, -1.62108055210833708e+09, 1.10684281682301447e+09, -4.95889784275030309e+08, 1.42062907797533095e+08, -2.44740627257387285e+07, 2.24376817792244943e+06, -8.40054336030240853e+04, 5.51335896122020586e+02, 8.14789096118312115e+08, -5.86648149205184723e+09, 1.86882075092958249e+10, -3.46320433881587779e+10, 4.12801855797539740e+10, -3.30265997498007231e+10, 1.79542137311556001e+10, -6.56329379261928433e+09, 1.55927986487925751e+09, -2.25105661889415278e+08, 1.73951075539781645e+07, -5.49842327572288687e+05, 3.03809051092238427e+03, -1.46792612476956167e+10, 1.14498237732025810e+11, -3.99096175224466498e+11, 8.19218669548577329e+11, -1.09837515608122331e+12, 1.00815810686538209e+12, -6.45364869245376503e+11, 2.87900649906150589e+11, -8.78670721780232657e+10, 1.76347306068349694e+10, -2.16716498322379509e+09, 1.43157876718888981e+08, -3.87183344257261262e+06, 1.82577554742931747e+04];
 alfa = [-4.44444444444444444e-03, -9.22077922077922078e-04, -8.84892884892884893e-05, 1.65927687832449737e-04, 2.46691372741792910e-04, 2.65995589346254780e-04, 2.61824297061500945e-04, 2.48730437344655609e-04, 2.32721040083232098e-04, 2.16362485712365082e-04, 2.00738858762752355e-04, 1.86267636637545172e-04, 1.73060775917876493e-04, 1.61091705929015752e-04, 1.50274774160908134e-04, 1.40503497391269794e-04, 1.31668816545922806e-04, 1.23667445598253261e-04, 1.16405271474737902e-04, 1.09798298372713369e-04, 1.03772410422992823e-04, 9.82626078369363448e-05, 9.32120517249503256e-05, 8.85710852478711718e-05, 8.42963105715700223e-05, 8.03497548407791151e-05, 7.66981345359207388e-05, 7.33122157481777809e-05, 7.01662625163141333e-05, 6.72375633790160292e-05, 6.93735541354588974e-04, 2.32241745182921654e-04, -1.41986273556691197e-05, -1.16444931672048640e-04, -1.50803558053048762e-04, -1.55121924918096223e-04, -1.46809756646465549e-04, -1.33815503867491367e-04, -1.19744975684254051e-04, -1.06184319207974020e-04, -9.37699549891194492e-05, -8.26923045588193274e-05, -7.29374348155221211e-05, -6.44042357721016283e-05, -5.69611566009369048e-05, -5.04731044303561628e-05, -4.48134868008882786e-05, -3.98688727717598864e-05, -3.55400532972042498e-05, -3.17414256609022480e-05, -2.83996793904174811e-05, -2.54522720634870566e-05, -2.28459297164724555e-05, -2.05352753106480604e-05, -1.84816217627666085e-05, -1.66519330021393806e-05, -1.50179412980119482e-05, -1.35554031379040526e-05, -1.22434746473858131e-05, -1.10641884811308169e-05, -3.54211971457743841e-04, -1.56161263945159416e-04, 3.04465503594936410e-05, 1.30198655773242693e-04, 1.67471106699712269e-04, 1.70222587683592569e-04, 1.56501427608594704e-04, 1.36339170977445120e-04, 1.14886692029825128e-04, 9.45869093034688111e-05, 7.64498419250898258e-05, 6.07570334965197354e-05, 4.74394299290508799e-05, 3.62757512005344297e-05, 2.69939714979224901e-05, 1.93210938247939253e-05, 1.30056674793963203e-05, 7.82620866744496661e-06, 3.59257485819351583e-06, 1.44040049814251817e-07, -2.65396769697939116e-06, -4.91346867098485910e-06, -6.72739296091248287e-06, -8.17269379678657923e-06, -9.31304715093561232e-06, -1.02011418798016441e-05, -1.08805962510592880e-05, -1.13875481509603555e-05, -1.17519675674556414e-05, -1.19987364870944141e-05, 3.78194199201772914e-04, 2.02471952761816167e-04, -6.37938506318862408e-05, -2.38598230603005903e-04, -3.10916256027361568e-04, -3.13680115247576316e-04, -2.78950273791323387e-04, -2.28564082619141374e-04, -1.75245280340846749e-04, -1.25544063060690348e-04, -8.22982872820208365e-05, -4.62860730588116458e-05, -1.72334302366962267e-05, 5.60690482304602267e-06, 2.31395443148286800e-05, 3.62642745856793957e-05, 4.58006124490188752e-05, 5.24595294959114050e-05, 5.68396208545815266e-05, 5.94349820393104052e-05, 6.06478527578421742e-05, 6.08023907788436497e-05, 6.01577894539460388e-05, 5.89199657344698500e-05, 5.72515823777593053e-05, 5.52804375585852577e-05, 5.31063773802880170e-05, 5.08069302012325706e-05, 4.84418647620094842e-05, 4.60568581607475370e-05, -6.91141397288294174e-04, -4.29976633058871912e-04, 1.83067735980039018e-04, 6.60088147542014144e-04, 8.75964969951185931e-04, 8.77335235958235514e-04, 7.49369585378990637e-04, 5.63832329756980918e-04, 3.68059319971443156e-04, 1.88464535514455599e-04, 3.70663057664904149e-05, -8.28520220232137023e-05, -1.72751952869172998e-04, -2.36314873605872983e-04, -2.77966150694906658e-04, -3.02079514155456919e-04, -3.12594712643820127e-04, -3.12872558758067163e-04, -3.05678038466324377e-04, -2.93226470614557331e-04, -2.77255655582934777e-04, -2.59103928467031709e-04, -2.39784014396480342e-04, -2.20048260045422848e-04, -2.00443911094971498e-04, -1.81358692210970687e-04, -1.63057674478657464e-04, -1.45712672175205844e-04, -1.29425421983924587e-04, -1.14245691942445952e-04, 1.92821964248775885e-03, 1.35592576302022234e-03, -7.17858090421302995e-04, -2.58084802575270346e-03, -3.49271130826168475e-03, -3.46986299340960628e-03, -2.82285233351310182e-03, -1.88103076404891354e-03, -8.89531718383947600e-04, 3.87912102631035228e-06, 7.28688540119691412e-04, 1.26566373053457758e-03, 1.62518158372674427e-03, 1.83203153216373172e-03, 1.91588388990527909e-03, 1.90588846755546138e-03, 1.82798982421825727e-03, 1.70389506421121530e-03, 1.55097127171097686e-03, 1.38261421852276159e-03, 1.20881424230064774e-03, 1.03676532638344962e-03, 8.71437918068619115e-04, 7.16080155297701002e-04, 5.72637002558129372e-04, 4.42089819465802277e-04, 3.24724948503090564e-04, 2.20342042730246599e-04, 1.28412898401353882e-04, 4.82005924552095464e-05];
 beta = [1.79988721413553309e-02, 5.59964911064388073e-03, 2.88501402231132779e-03, 1.80096606761053941e-03, 1.24753110589199202e-03, 9.22878876572938311e-04, 7.14430421727287357e-04, 5.71787281789704872e-04, 4.69431007606481533e-04, 3.93232835462916638e-04, 3.34818889318297664e-04, 2.88952148495751517e-04, 2.52211615549573284e-04, 2.22280580798883327e-04, 1.97541838033062524e-04, 1.76836855019718004e-04, 1.59316899661821081e-04, 1.44347930197333986e-04, 1.31448068119965379e-04, 1.20245444949302884e-04, 1.10449144504599392e-04, 1.01828770740567258e-04, 9.41998224204237509e-05, 8.74130545753834437e-05, 8.13466262162801467e-05, 7.59002269646219339e-05, 7.09906300634153481e-05, 6.65482874842468183e-05, 6.25146958969275078e-05, 5.88403394426251749e-05, -1.49282953213429172e-03, -8.78204709546389328e-04, -5.02916549572034614e-04, -2.94822138512746025e-04, -1.75463996970782828e-04, -1.04008550460816434e-04, -5.96141953046457895e-05, -3.12038929076098340e-05, -1.26089735980230047e-05, -2.42892608575730389e-07, 8.05996165414273571e-06, 1.36507009262147391e-05, 1.73964125472926261e-05, 1.98672978842133780e-05, 2.14463263790822639e-05, 2.23954659232456514e-05, 2.28967783814712629e-05, 2.30785389811177817e-05, 2.30321976080909144e-05, 2.28236073720348722e-05, 2.25005881105292418e-05, 2.20981015361991429e-05, 2.16418427448103905e-05, 2.11507649256220843e-05, 2.06388749782170737e-05, 2.01165241997081666e-05, 1.95913450141179244e-05, 1.90689367910436740e-05, 1.85533719641636667e-05, 1.80475722259674218e-05, 5.52213076721292790e-04, 4.47932581552384646e-04, 2.79520653992020589e-04, 1.52468156198446602e-04, 6.93271105657043598e-05, 1.76258683069991397e-05, -1.35744996343269136e-05, -3.17972413350427135e-05, -4.18861861696693365e-05, -4.69004889379141029e-05, -4.87665447413787352e-05, -4.87010031186735069e-05, -4.74755620890086638e-05, -4.55813058138628452e-05, -4.33309644511266036e-05, -4.09230193157750364e-05, -3.84822638603221274e-05, -3.60857167535410501e-05, -3.37793306123367417e-05, -3.15888560772109621e-05, -2.95269561750807315e-05, -2.75978914828335759e-05, -2.58006174666883713e-05, -2.41308356761280200e-05, -2.25823509518346033e-05, -2.11479656768912971e-05, -1.98200638885294927e-05, -1.85909870801065077e-05, -1.74532699844210224e-05, -1.63997823854497997e-05, -4.74617796559959808e-04, -4.77864567147321487e-04, -3.20390228067037603e-04, -1.61105016119962282e-04, -4.25778101285435204e-05, 3.44571294294967503e-05, 7.97092684075674924e-05, 1.03138236708272200e-04, 1.12466775262204158e-04, 1.13103642108481389e-04, 1.08651634848774268e-04, 1.01437951597661973e-04, 9.29298396593363896e-05, 8.40293133016089978e-05, 7.52727991349134062e-05, 6.69632521975730872e-05, 5.92564547323194704e-05, 5.22169308826975567e-05, 4.58539485165360646e-05, 4.01445513891486808e-05, 3.50481730031328081e-05, 3.05157995034346659e-05, 2.64956119950516039e-05, 2.29363633690998152e-05, 1.97893056664021636e-05, 1.70091984636412623e-05, 1.45547428261524004e-05, 1.23886640995878413e-05, 1.04775876076583236e-05, 8.79179954978479373e-06, 7.36465810572578444e-04, 8.72790805146193976e-04, 6.22614862573135066e-04, 2.85998154194304147e-04, 3.84737672879366102e-06, -1.87906003636971558e-04, -2.97603646594554535e-04, -3.45998126832656348e-04, -3.53382470916037712e-04, -3.35715635775048757e-04, -3.04321124789039809e-04, -2.66722723047612821e-04, -2.27654214122819527e-04, -1.89922611854562356e-04, -1.55058918599093870e-04, -1.23778240761873630e-04, -9.62926147717644187e-05, -7.25178327714425337e-05, -5.22070028895633801e-05, -3.50347750511900522e-05, -2.06489761035551757e-05, -8.70106096849767054e-06, 1.13698686675100290e-06, 9.16426474122778849e-06, 1.56477785428872620e-05, 2.08223629482466847e-05, 2.48923381004595156e-05, 2.80340509574146325e-05, 3.03987774629861915e-05, 3.21156731406700616e-05, -1.80182191963885708e-03, -2.43402962938042533e-03, -1.83422663549856802e-03, -7.62204596354009765e-04, 2.39079475256927218e-04, 9.49266117176881141e-04, 1.34467449701540359e-03, 1.48457495259449178e-03, 1.44732339830617591e-03, 1.30268261285657186e-03, 1.10351597375642682e-03, 8.86047440419791759e-04, 6.73073208165665473e-04, 4.77603872856582378e-04, 3.05991926358789362e-04, 1.60315694594721630e-04, 4.00749555270613286e-05, -5.66607461635251611e-05, -1.32506186772982638e-04, -1.90296187989614057e-04, -2.32811450376937408e-04, -2.62628811464668841e-04, -2.82050469867598672e-04, -2.93081563192861167e-04, -2.97435962176316616e-04, -2.96557334239348078e-04, -2.91647363312090861e-04, -2.83696203837734166e-04, -2.73512317095673346e-04, -2.61750155806768580e-04, 6.38585891212050914e-03, 9.62374215806377941e-03, 7.61878061207001043e-03, 2.83219055545628054e-03, -2.09841352012720090e-03, -5.73826764216626498e-03, -7.70804244495414620e-03, -8.21011692264844401e-03, -7.65824520346905413e-03, -6.47209729391045177e-03, -4.99132412004966473e-03, -3.45612289713133280e-03, -2.01785580014170775e-03, -7.59430686781961401e-04, 2.84173631523859138e-04, 1.10891667586337403e-03, 1.72901493872728771e-03, 2.16812590802684701e-03, 2.45357710494539735e-03, 2.61281821058334862e-03, 2.67141039656276912e-03, 2.65203073395980430e-03, 2.57411652877287315e-03, 2.45389126236094427e-03, 2.30460058071795494e-03, 2.13684837686712662e-03, 1.95896528478870911e-03, 1.77737008679454412e-03, 1.59690280765839059e-03, 1.42111975664438546e-03];
 gama = [6.29960524947436582e-01, 2.51984209978974633e-01, 1.54790300415655846e-01, 1.10713062416159013e-01, 8.57309395527394825e-02, 6.97161316958684292e-02, 5.86085671893713576e-02, 5.04698873536310685e-02, 4.42600580689154809e-02, 3.93720661543509966e-02, 3.54283195924455368e-02, 3.21818857502098231e-02, 2.94646240791157679e-02, 2.71581677112934479e-02, 2.51768272973861779e-02, 2.34570755306078891e-02, 2.19508390134907203e-02, 2.06210828235646240e-02, 1.94388240897880846e-02, 1.83810633800683158e-02, 1.74293213231963172e-02, 1.65685837786612353e-02, 1.57865285987918445e-02, 1.50729501494095594e-02, 1.44193250839954639e-02, 1.38184805735341786e-02, 1.32643378994276568e-02, 1.27517121970498651e-02, 1.22761545318762767e-02, 1.18338262398482403e-02];
 ex1 = 3.33333333333333333e-01;
 ex2 = 6.66666666666666667e-01;
 hpi = 1.57079632679489662;
 gpi = 3.14159265358979324;
 thpi = 4.71238898038468986;
 zeror = 0;
 zeroi = 0;
 coner = 1;
 conei = 0;


 var goToLabel = 0;
 mainExecutionLoop: while (true) {
 switch (goToLabel) {
 case 0:
 rfnu = 1.0 / fnu;
 // c-----------------------------------------------------------------------
 // c overflow test (z/fnu too small)
 // c-----------------------------------------------------------------------
 test = (0, _d1mach.d1mach)(1) * 1.0e+3;
 ac = fnu * test;
 if (Math.abs(zr) > ac || Math.abs(zi) > ac) {
 goToLabel = 15;break;
 }
 zeta1r = 2.0 * Math.abs(Math.log(test)) + fnu;
 zeta1i = 0.0;
 zeta2r = fnu;
 zeta2i = 0.0;
 phir = 1.0;
 phii = 0.0;
 argr = 1.0;
 argi = 0.0;
 break mainExecutionLoop;
 case 15:
 zbr = zr * rfnu;
 zbi = zi * rfnu;
 rfnu2 = rfnu * rfnu;
 // c-----------------------------------------------------------------------
 // c compute in the fourth quadrant
 // c-----------------------------------------------------------------------
 fn13 = fnu ** ex1;
 fn23 = fn13 * fn13;
 rfn13 = 1.0 / fn13;
 w2r = coner - zbr * zbr + zbi * zbi;
 w2i = conei - zbr * zbi - zbr * zbi;
 aw2 = (0, _zabs.azabs)(w2r, w2i);
 if (aw2 > 0.25) {
 goToLabel = 130;break;
 }
 // c-----------------------------------------------------------------------
 // c power series for cabs(w2) <= 0.25
 // c-----------------------------------------------------------------------
 k = 1;
 pr[0] = coner;
 pi[0] = conei;
 sumar = gama[0];
 sumai = zeroi;
 ap[0] = 1.0;
 if (aw2 < tol) {
 goToLabel = 20;break;
 }
 // do 10 k=2,30
 for (k = 2; k <= 30; k++) {
 pr[k - 1] = pr[k - 2] * w2r - pi[k - 2] * w2i;
 pi[k - 1] = pr[k - 2] * w2i + pi[k - 2] * w2r;
 sumar = sumar + pr[k - 1] * gama[k - 1];
 sumai = sumai + pi[k - 1] * gama[k - 1];
 ap[k - 1] = ap[k - 2] * aw2;
 if (ap[k - 1] < tol) {
 goToLabel = 20;break;
 }
 }
 // 10 continue
 if (goToLabel < 20) k = 30;
 case 20:
 kmax = k;
 zetar = w2r * sumar - w2i * sumai;
 zetai = w2r * sumai + w2i * sumar;
 argr = zetar * fn23;
 argi = zetai * fn23;

 var _azsqrt = (0, _zsqrt.azsqrt)(sumar, sumai);

 var _azsqrt2 = _slicedToArray(_azsqrt, 2);

 zar = _azsqrt2[0];
 zai = _azsqrt2[1];

 var _azsqrt3 = (0, _zsqrt.azsqrt)(w2r, w2i);

 var _azsqrt4 = _slicedToArray(_azsqrt3, 2);

 str = _azsqrt4[0];
 sti = _azsqrt4[1];

 zeta2r = str * fnu;
 zeta2i = sti * fnu;
 str = coner + ex2 * (zetar * zar - zetai * zai);
 sti = conei + ex2 * (zetar * zai + zetai * zar);
 zeta1r = str * zeta2r - sti * zeta2i;
 zeta1i = str * zeta2i + sti * zeta2r;
 zar = zar + zar;
 zai = zai + zai;

 var _azsqrt5 = (0, _zsqrt.azsqrt)(zar, zai);

 var _azsqrt6 = _slicedToArray(_azsqrt5, 2);

 str = _azsqrt6[0];
 sti = _azsqrt6[1];

 phir = str * rfn13;
 phii = sti * rfn13;
 if (ipmtr === 1) {
 goToLabel = 120;break;
 }
 // c-----------------------------------------------------------------------
 // c sum series for asum and bsum
 // c-----------------------------------------------------------------------
 sumbr = zeror;
 sumbi = zeroi;
 // do 30 k=1,kmax
 for (k = 1; k <= kmax; k++) {
 sumbr = sumbr + pr[k - 1] * beta[k - 1];
 sumbi = sumbi + pi[k - 1] * beta[k - 1];
 }
 // 30 continue
 asumr = zeror;
 asumi = zeroi;
 bsumr = sumbr;
 bsumi = sumbi;
 l1 = 0;
 l2 = 30;
 btol = tol * (Math.abs(bsumr) + Math.abs(bsumi));
 atol = tol;
 pp = 1.0;
 ias = 0;
 ibs = 0;
 if (rfnu2 < tol) {
 goToLabel = 110;break;
 }
 // do 100 is=2,7
 for (is = 2; is <= 7; is++) {
 atol = atol / rfnu2;
 pp = pp * rfnu2;
 if (ias === 1) {
 // go to 60
 } else {
 sumar = zeror;
 sumai = zeroi;
 // do 40 k=1,kmax
 for (k = 1; k <= kmax; k++) {
 m = l1 + k;
 sumar = sumar + pr[k - 1] * alfa[m - 1];
 sumai = sumai + pi[k - 1] * alfa[m - 1];
 if (ap[k - 1] < atol) break;
 }
 // 40 continue
 // 50 continue
 asumr = asumr + sumar * pp;
 asumi = asumi + sumai * pp;
 if (pp < tol) ias = 1;
 }
 // 60 continue
 if (ibs === 1) {
 // go to 90
 } else {
 sumbr = zeror;
 sumbi = zeroi;
 // do 70 k=1,kmax
 for (k = 1; k <= kmax; k++) {
 m = l2 + k;
 sumbr = sumbr + pr[k - 1] * beta[m - 1];
 sumbi = sumbi + pi[k - 1] * beta[m - 1];
 if (ap[k - 1] < atol) break;
 }
 // 70 continue
 // 80 continue
 bsumr = bsumr + sumbr * pp;
 bsumi = bsumi + sumbi * pp;
 if (pp < btol) ibs = 1;
 }
 // 90 continue
 if (ias === 1 && ibs === 1) break;
 l1 = l1 + 30;
 l2 = l2 + 30;
 }
 // 100 continue
 case 110:
 asumr = asumr + coner;
 pp = rfnu * rfn13;
 bsumr = bsumr * pp;
 bsumi = bsumi * pp;
 case 120:
 break mainExecutionLoop;
 // c-----------------------------------------------------------------------
 // c cabs(w2) > 0.25
 // c-----------------------------------------------------------------------
 case 130:
 var _azsqrt7 = (0, _zsqrt.azsqrt)(w2r, w2i);

 var _azsqrt8 = _slicedToArray(_azsqrt7, 2);

 wr = _azsqrt8[0];
 wi = _azsqrt8[1];

 if (wr < 0.0) wr = 0.0;
 if (wi < 0.0) wi = 0.0;
 str = coner + wr;
 sti = wi;

 var _zdiv = (0, _zdiv9.zdiv)(str, sti, zbr, zbi);

 var _zdiv2 = _slicedToArray(_zdiv, 2);

 zar = _zdiv2[0];
 zai = _zdiv2[1];

 var _azlog = (0, _zlog.azlog)(zar, zai);

 var _azlog2 = _slicedToArray(_azlog, 2);

 zcr = _azlog2[0];
 zci = _azlog2[1];

 if (zci < 0.0) zci = 0.0;
 if (zci > hpi) zci = hpi;
 if (zcr < 0.0) zcr = 0.0;
 zthr = (zcr - wr) * 1.5;
 zthi = (zci - wi) * 1.5;
 zeta1r = zcr * fnu;
 zeta1i = zci * fnu;
 zeta2r = wr * fnu;
 zeta2i = wi * fnu;
 azth = (0, _zabs.azabs)(zthr, zthi);
 ang = thpi;
 if (zthr >= 0.0 && zthi < 0.0) {
 goToLabel = 140;break;
 }
 ang = hpi;
 if (zthr === 0.0) {
 goToLabel = 140;break;
 }
 ang = Math.atan(zthi / zthr);
 if (zthr < 0.0) ang = ang + gpi;
 case 140:
 pp = azth ** ex2;
 ang = ang * ex2;
 zetar = pp * Math.cos(ang);
 zetai = pp * Math.sin(ang);
 if (zetai < 0.0) zetai = 0.0;
 argr = zetar * fn23;
 argi = zetai * fn23;

 var _zdiv3 = (0, _zdiv9.zdiv)(zthr, zthi, zetar, zetai);

 var _zdiv4 = _slicedToArray(_zdiv3, 2);

 rtztr = _zdiv4[0];
 rtzti = _zdiv4[1];

 var _zdiv5 = (0, _zdiv9.zdiv)(rtztr, rtzti, wr, wi);

 var _zdiv6 = _slicedToArray(_zdiv5, 2);

 zar = _zdiv6[0];
 zai = _zdiv6[1];

 tzar = zar + zar;
 tzai = zai + zai;

 var _azsqrt9 = (0, _zsqrt.azsqrt)(tzar, tzai);

 var _azsqrt10 = _slicedToArray(_azsqrt9, 2);

 str = _azsqrt10[0];
 sti = _azsqrt10[1];

 phir = str * rfn13;
 phii = sti * rfn13;
 if (ipmtr === 1) {
 goToLabel = 120;break;
 }
 raw = 1.0 / Math.sqrt(aw2);
 str = wr * raw;
 sti = -wi * raw;
 tfnr = str * rfnu * raw;
 tfni = sti * rfnu * raw;
 razth = 1.0 / azth;
 str = zthr * razth;
 sti = -zthi * razth;
 rzthr = str * razth * rfnu;
 rzthi = sti * razth * rfnu;
 zcr = rzthr * ar[1];
 zci = rzthi * ar[1];
 raw2 = 1.0 / aw2;
 str = w2r * raw2;
 sti = -w2i * raw2;
 t2r = str * raw2;
 t2i = sti * raw2;
 str = t2r * c[1] + c[2];
 sti = t2i * c[1];
 upr[1] = str * tfnr - sti * tfni;
 upi[1] = str * tfni + sti * tfnr;
 bsumr = upr[1] + zcr;
 bsumi = upi[1] + zci;
 asumr = zeror;
 asumi = zeroi;
 if (rfnu < tol) {
 goToLabel = 220;break;
 }
 przthr = rzthr;
 przthi = rzthi;
 ptfnr = tfnr;
 ptfni = tfni;
 upr[0] = coner;
 upi[0] = conei;
 pp = 1.0;
 btol = tol * (Math.abs(bsumr) + Math.abs(bsumi));
 ks = 0;
 kp1 = 2;
 l = 3;
 ias = 0;
 ibs = 0;
 // do 210 lr=2,12,2
 for (lr = 2; lr <= 12; lr += 2) {
 lrp1 = lr + 1;
 // c-----------------------------------------------------------------------
 // c compute two additional cr, dr, and up for two more terms in
 // c next suma and sumb
 // c-----------------------------------------------------------------------
 // do 160 k=lr,lrp1
 for (k = lr; k <= lrp1; k++) {
 ks = ks + 1;
 kp1 = kp1 + 1;
 l = l + 1;
 zar = c[l - 1];
 zai = zeroi;
 // do 150 j=2,kp1
 for (j = 2; j <= kp1; j++) {
 l = l + 1;
 str = zar * t2r - t2i * zai + c[l - 1];
 zai = zar * t2i + zai * t2r;
 zar = str;
 }
 // 150 continue
 str = ptfnr * tfnr - ptfni * tfni;
 ptfni = ptfnr * tfni + ptfni * tfnr;
 ptfnr = str;
 upr[kp1 - 1] = ptfnr * zar - ptfni * zai;
 upi[kp1 - 1] = ptfni * zar + ptfnr * zai;
 crr[ks - 1] = przthr * br[ks];
 cri[ks - 1] = przthi * br[ks];
 str = przthr * rzthr - przthi * rzthi;
 przthi = przthr * rzthi + przthi * rzthr;
 przthr = str;
 drr[ks - 1] = przthr * ar[ks + 1];
 dri[ks - 1] = przthi * ar[ks + 1];
 }
 pp = pp * rfnu2;
 if (ias === 1) {
 // go to 180
 } else {
 sumar = upr[lrp1 - 1];
 sumai = upi[lrp1 - 1];
 ju = lrp1;
 // do 170 jr=1,lr
 for (jr = 1; jr <= lr; jr++) {
 ju = ju - 1;
 sumar = sumar + crr[jr - 1] * upr[ju - 1] - cri[jr - 1] * upi[ju - 1];
 sumai = sumai + crr[jr - 1] * upi[ju - 1] + cri[jr - 1] * upr[ju - 1];
 }
 // 170 continue
 asumr = asumr + sumar;
 asumi = asumi + sumai;
 test = Math.abs(sumar) + Math.abs(sumai);
 if (pp < tol && test < tol) ias = 1;
 }
 // 180 continue
 if (ibs === 1) {
 // go to 200
 } else {
 sumbr = upr[lr + 1] + upr[lrp1 - 1] * zcr - upi[lrp1 - 1] * zci;
 sumbi = upi[lr + 1] + upr[lrp1 - 1] * zci + upi[lrp1 - 1] * zcr;
 ju = lrp1;
 // do 190 jr=1,lr
 for (jr = 1; jr <= lr; jr++) {
 ju = ju - 1;
 sumbr = sumbr + drr[jr - 1] * upr[ju - 1] - dri[jr - 1] * upi[ju - 1];
 sumbi = sumbi + drr[jr - 1] * upi[ju - 1] + dri[jr - 1] * upr[ju - 1];
 }
 // 190 continue
 bsumr = bsumr + sumbr;
 bsumi = bsumi + sumbi;
 test = Math.abs(sumbr) + Math.abs(sumbi);
 if (pp < btol && test < btol) ibs = 1;
 }
 // 200 continue
 if (ias === 1 && ibs === 1) break;
 }
 // 210 continue
 case 220:
 asumr = asumr + coner;
 str = -bsumr * rfn13;
 sti = -bsumi * rfn13;

 var _zdiv7 = (0, _zdiv9.zdiv)(str, sti, rtztr, rtzti);

 var _zdiv8 = _slicedToArray(_zdiv7, 2);

 bsumr = _zdiv8[0];
 bsumi = _zdiv8[1];

 goToLabel = 120;break;
 default:
 break mainExecutionLoop;
 }
 }

 return [phir, phii, argr, argi, zeta1r, zeta1i, zeta2r, zeta2i, asumr, asumi, bsumr, bsumi];
}
},{"../../utils/fortran-utils/d1mach.js":91,"./zabs.js":11,"./zdiv.js":26,"./zlog.js":29,"./zsqrt.js":36}],39:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});

var _slicedToArray = function () { function sliceIterator(arr, i) { var _arr = []; var _n = true; var _d = false; var _e = undefined; try { for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) { _arr.push(_s.value); if (i && _arr.length === i) break; } } catch (err) { _d = true; _e = err; } finally { try { if (!_n && _i["return"]) _i["return"](); } finally { if (_d) throw _e; } } return _arr; } return function (arr, i) { if (Array.isArray(arr)) { return arr; } else if (Symbol.iterator in Object(arr)) { return sliceIterator(arr, i); } else { throw new TypeError("Invalid attempt to destructure non-iterable instance"); } }; }(); /* eslint-disable no-labels */
/* eslint-disable no-fallthrough */
// SUBROUTINE ZUNI1(ZR, ZI, FNU, KODE, N, YR, YI, NZ, NLAST, FNUL,
// * TOL, ELIM, ALIM)
// ***BEGIN PROLOGUE ZUNI1
// ***REFER TO ZBESI,ZBESK
//
// ZUNI1 COMPUTES I(FNU,Z) BY MEANS OF THE UNIFORM ASYMPTOTIC
// EXPANSION FOR I(FNU,Z) IN -PI/3.LE.ARG Z.LE.PI/3.
//
// FNUL IS THE SMALLEST ORDER PERMITTED FOR THE ASYMPTOTIC
// EXPANSION. NLAST=0 MEANS ALL OF THE Y VALUES WERE SET.
// NLAST.NE.0 IS THE NUMBER LEFT TO BE COMPUTED BY ANOTHER
// FORMULA FOR ORDERS FNU TO FNU+NLAST-1 BECAUSE FNU+NLAST-1.LT.FNUL.
// Y(I)=CZERO FOR I=NLAST+1,N
//
// ***ROUTINES CALLED ZUCHK,ZUNIK,ZUOIK,D1MACH,AZABS
// ***END PROLOGUE ZUNI1


exports.zuni1 = zuni1;

var _zuchk = require('./zuchk.js');

var _zunik5 = require('./zunik.js');

var _zuoik = require('./zuoik.js');

var _d1mach = require('../../utils/fortran-utils/d1mach.js');

var _zabs = require('./zabs.js');

function zuni1(zr, zi, fnu, kode, n, yr, yi, fnul, tol, elim, alim) {
 var aphi = void 0,
 ascle = void 0,
 bry = void 0,
 coner = void 0,
 crsc = void 0,
 cscl = void 0,
 csrr = void 0,
 cssr = void 0,
 c1r = void 0,
 c2i = void 0,
 c2m = void 0,
 c2r = void 0,
 fn = void 0,
 phii = void 0,
 phir = void 0,
 rast = void 0,
 rs1 = void 0,
 rzi = void 0,
 rzr = void 0,
 sti = void 0,
 str = void 0,
 sumi = void 0,
 sumr = void 0,
 s1i = void 0,
 s1r = void 0,
 s2i = void 0,
 s2r = void 0,
 zeroi = void 0,
 zeror = void 0,
 zeta1i = void 0,
 zeta1r = void 0,
 zeta2i = void 0,
 zeta2r = void 0,
 cyr = void 0,
 cyi = void 0,
 i = void 0,
 iflag = void 0,
 init = void 0,
 k = void 0,
 m = void 0,
 nd = void 0,
 nlast = void 0,
 nn = void 0,
 nuf = void 0,
 nw = void 0,
 nz = void 0;
 bry = new Array(3);
 cssr = new Array(3);
 csrr = new Array(3);
 cyr = new Array(2);
 cyi = new Array(2);
 var goToLabel = 0;
 mainExecutionLoop: while (true) {
 switch (goToLabel) {
 case 0:
 zeror = 0.0;
 zeroi = 0.0;
 coner = 1.0;

 nz = 0;
 nd = n;
 nlast = 0;
 // c-----------------------------------------------------------------------
 // c computed values with exponents between alim and elim in mag-
 // c nitude are scaled to keep intermediate arithmetic on scale,
 // c exp(alim)=exp(elim)*tol
 // c-----------------------------------------------------------------------
 cscl = 1.0 / tol;
 crsc = tol;
 cssr[0] = cscl;
 cssr[1] = coner;
 cssr[2] = crsc;
 csrr[0] = crsc;
 csrr[1] = coner;
 csrr[2] = cscl;
 bry[0] = 1.0e+3 * (0, _d1mach.d1mach)(1) / tol;
 // c-----------------------------------------------------------------------
 // c check for underflow and overflow on first member
 // c-----------------------------------------------------------------------
 fn = Math.max(fnu, 1.0);
 init = 0;

 var _zunik = (0, _zunik5.zunik)(zr, zi, fn, 1, 1, tol, init);

 var _zunik2 = _slicedToArray(_zunik, 8);

 phir = _zunik2[0];
 phii = _zunik2[1];
 zeta1r = _zunik2[2];
 zeta1i = _zunik2[3];
 zeta2r = _zunik2[4];
 zeta2i = _zunik2[5];
 sumr = _zunik2[6];
 sumi = _zunik2[7];

 if (kode === 1) {
 goToLabel = 10;break;
 }
 str = zr + zeta2r;
 sti = zi + zeta2i;
 rast = fn / (0, _zabs.azabs)(str, sti);
 str = str * rast * rast;
 sti = -sti * rast * rast;
 s1r = -zeta1r + str;
 s1i = -zeta1i + sti;
 goToLabel = 20;break;
 case 10:
 s1r = -zeta1r + zeta2r;
 s1i = -zeta1i + zeta2i;
 case 20:
 rs1 = s1r;
 if (Math.abs(rs1) > elim) {
 goToLabel = 130;break;
 }
 case 30:
 nn = Math.min(2, nd);
 // do 80 i=1,nn
 forLoop80: for (i = 1; i <= nn; i++) {
 fn = fnu + (nd - i);
 init = 0;

 var _zunik3 = (0, _zunik5.zunik)(zr, zi, fn, 1, 0, tol, init);

 var _zunik4 = _slicedToArray(_zunik3, 8);

 phir = _zunik4[0];
 phii = _zunik4[1];
 zeta1r = _zunik4[2];
 zeta1i = _zunik4[3];
 zeta2r = _zunik4[4];
 zeta2i = _zunik4[5];
 sumr = _zunik4[6];
 sumi = _zunik4[7];

 if (kode === 1) {
 s1r = -zeta1r + zeta2r;
 s1i = -zeta1i + zeta2i;
 } else {
 str = zr + zeta2r;
 sti = zi + zeta2i;
 rast = fn / (0, _zabs.azabs)(str, sti);
 str = str * rast * rast;
 sti = -sti * rast * rast;
 s1r = -zeta1r + str;
 s1i = -zeta1i + sti + zi;
 }
 // c-----------------------------------------------------------------------
 // c test for underflow and overflow
 // c-----------------------------------------------------------------------
 rs1 = s1r;
 if (Math.abs(rs1) > elim) {
 goToLabel = 110;break forLoop80;
 }
 if (i === 1) iflag = 2;
 if (Math.abs(rs1) < alim) {
 // go to 60
 } else {
 // c-----------------------------------------------------------------------
 // c refine test and scale
 // c-----------------------------------------------------------------------
 aphi = (0, _zabs.azabs)(phir, phii);
 rs1 = rs1 + Math.log(aphi);
 if (Math.abs(rs1) > elim) {
 goToLabel = 110;break forLoop80;
 }
 if (i === 1) iflag = 1;
 if (rs1 < 0.0) {
 // go to 60
 } else {
 if (i === 1) iflag = 3;
 }
 }
 // 60 continue
 // c-----------------------------------------------------------------------
 // c scale s1 if cabs(s1) < ascle
 // c-----------------------------------------------------------------------
 s2r = phir * sumr - phii * sumi;
 s2i = phir * sumi + phii * sumr;
 str = Math.exp(s1r) * cssr[iflag - 1];
 s1r = str * Math.cos(s1i);
 s1i = str * Math.sin(s1i);
 str = s2r * s1r - s2i * s1i;
 s2i = s2r * s1i + s2i * s1r;
 s2r = str;
 if (iflag !== 1) {
 // go to 70
 } else {
 nw = (0, _zuchk.zuchk)(s2r, s2i, bry[0], tol);
 if (nw !== 0) {
 goToLabel = 110;break forLoop80;
 }
 }
 // 70 continue
 cyr[i - 1] = s2r;
 cyi[i - 1] = s2i;
 m = nd - i + 1;
 yr[m - 1] = s2r * csrr[iflag - 1];
 yi[m - 1] = s2i * csrr[iflag - 1];
 }
 // 80 continue
 if (goToLabel > 80) {
 break;
 }
 if (nd <= 2) {
 goToLabel = 100;break;
 }
 rast = 1.0 / (0, _zabs.azabs)(zr, zi);
 str = zr * rast;
 sti = -zi * rast;
 rzr = (str + str) * rast;
 rzi = (sti + sti) * rast;
 bry[1] = 1.0 / bry[0];
 bry[2] = (0, _d1mach.d1mach)(2);
 s1r = cyr[0];
 s1i = cyi[0];
 s2r = cyr[1];
 s2i = cyi[1];
 c1r = csrr[iflag - 1];
 ascle = bry[iflag - 1];
 k = nd - 2;
 fn = k;

 // do 90 i=3,nd
 for (i = 3; i <= nd; i++) {
 c2r = s2r;
 c2i = s2i;
 s2r = s1r + (fnu + fn) * (rzr * c2r - rzi * c2i);
 s2i = s1i + (fnu + fn) * (rzr * c2i + rzi * c2r);
 s1r = c2r;
 s1i = c2i;
 c2r = s2r * c1r;
 c2i = s2i * c1r;
 yr[k - 1] = c2r;
 yi[k - 1] = c2i;
 k = k - 1;
 fn = fn - 1.0;
 if (iflag >= 3) break;
 str = Math.abs(c2r);
 sti = Math.abs(c2i);
 c2m = Math.max(str, sti);
 if (c2m <= ascle) break;
 iflag = iflag + 1;
 ascle = bry[iflag - 1];
 s1r = s1r * c1r;
 s1i = s1i * c1r;
 s2r = c2r;
 s2i = c2i;
 s1r = s1r * cssr[iflag - 1];
 s1i = s1i * cssr[iflag - 1];
 s2r = s2r * cssr[iflag - 1];
 s2i = s2i * cssr[iflag - 1];
 c1r = csrr[iflag - 1];
 }
 // 90 continue
 case 100:
 break mainExecutionLoop;
 case 110:
 // c-----------------------------------------------------------------------
 // c set underflow and update parameters
 // c-----------------------------------------------------------------------
 if (rs1 > 0.0) {
 goToLabel = 120;break;
 }
 yr[nd - 1] = zeror;
 yi[nd - 1] = zeroi;
 nz = nz + 1;
 nd = nd - 1;
 if (nd === 0) {
 goToLabel = 100;break;
 }
 nuf = (0, _zuoik.zuoik)(zr, zi, fnu, kode, 1, nd, yr, yi, tol, elim, alim);
 if (nuf < 0) {
 goToLabel = 120;break;
 }
 nd = nd - nuf;
 nz = nz + nuf;
 if (nd === 0) {
 goToLabel = 100;break;
 }
 fn = fnu + (nd - 1);
 if (fn >= fnul) {
 goToLabel = 30;break;
 }
 nlast = nd;
 break mainExecutionLoop;
 case 120:
 nz = -1;
 break mainExecutionLoop;
 case 130:
 if (rs1 > 0.0) {
 goToLabel = 120;break;
 }
 nz = n;
 // do 140 i=1,n
 for (i = 1; i <= n; i++) {
 yr[i - 1] = zeror;
 yi[i - 1] = zeroi;
 }
 // 140 continue
 default:
 break mainExecutionLoop;
 }
 }

 return [nz, nlast];
}
},{"../../utils/fortran-utils/d1mach.js":91,"./zabs.js":11,"./zuchk.js":37,"./zunik.js":41,"./zuoik.js":44}],40:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});

var _slicedToArray = function () { function sliceIterator(arr, i) { var _arr = []; var _n = true; var _d = false; var _e = undefined; try { for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) { _arr.push(_s.value); if (i && _arr.length === i) break; } } catch (err) { _d = true; _e = err; } finally { try { if (!_n && _i["return"]) _i["return"](); } finally { if (_d) throw _e; } } return _arr; } return function (arr, i) { if (Array.isArray(arr)) { return arr; } else if (Symbol.iterator in Object(arr)) { return sliceIterator(arr, i); } else { throw new TypeError("Invalid attempt to destructure non-iterable instance"); } }; }(); /* eslint-disable no-labels */
/* eslint-disable no-fallthrough */
// SUBROUTINE ZUNI2(ZR, ZI, FNU, KODE, N, YR, YI, NZ, NLAST, FNUL,
// * TOL, ELIM, ALIM)
// ***BEGIN PROLOGUE ZUNI2
// ***REFER TO ZBESI,ZBESK
//
// ZUNI2 COMPUTES I(FNU,Z) IN THE RIGHT HALF PLANE BY MEANS OF
// UNIFORM ASYMPTOTIC EXPANSION FOR J(FNU,ZN) WHERE ZN IS Z*I
// OR -Z*I AND ZN IS IN THE RIGHT HALF PLANE ALSO.
//
// FNUL IS THE SMALLEST ORDER PERMITTED FOR THE ASYMPTOTIC
// EXPANSION. NLAST=0 MEANS ALL OF THE Y VALUES WERE SET.
// NLAST.NE.0 IS THE NUMBER LEFT TO BE COMPUTED BY ANOTHER
// FORMULA FOR ORDERS FNU TO FNU+NLAST-1 BECAUSE FNU+NLAST-1.LT.FNUL.
// Y(I)=CZERO FOR I=NLAST+1,N
//
// ***ROUTINES CALLED ZAIRY,ZUCHK,ZUNHJ,ZUOIK,D1MACH,AZABS
// ***END PROLOGUE ZUNI2


exports.zuni2 = zuni2;

var _zairy5 = require('./zairy.js');

var _zuchk = require('./zuchk.js');

var _zunhj5 = require('./zunhj.js');

var _zuoik = require('./zuoik.js');

var _zabs = require('./zabs.js');

var _d1mach = require('../../utils/fortran-utils/d1mach.js');

function zuni2(zr, zi, fnu, kode, n, yr, yi, fnul, tol, elim, alim) {
 var aarg = void 0,
 aic = void 0,
 aii = void 0,
 air = void 0,
 ang = void 0,
 aphi = void 0,
 argi = void 0,
 argr = void 0,
 ascle = void 0,
 asumi = void 0,
 asumr = void 0,
 bry = void 0,
 bsumi = void 0,
 bsumr = void 0,
 cidi = void 0,
 cipi = void 0,
 cipr = void 0,
 coner = void 0,
 crsc = void 0,
 cscl = void 0,
 csrr = void 0,
 cssr = void 0,
 c1r = void 0,
 c2i = void 0,
 c2m = void 0,
 c2r = void 0,
 daii = void 0,
 dair = void 0,
 fn = void 0,
 hpi = void 0,
 phii = void 0,
 phir = void 0,
 rast = void 0,
 raz = void 0,
 rs1 = void 0,
 rzi = void 0,
 rzr = void 0,
 sti = void 0,
 str = void 0,
 s1i = void 0,
 s1r = void 0,
 s2i = void 0,
 s2r = void 0,
 zbi = void 0,
 zbr = void 0,
 zeroi = void 0,
 zeror = void 0,
 zeta1i = void 0,
 zeta1r = void 0,
 zeta2i = void 0,
 zeta2r = void 0,
 zni = void 0,
 znr = void 0,
 cyr = void 0,
 cyi = void 0,
 car = void 0,
 sar = void 0,
 i = void 0,
 iflag = void 0,
 index = void 0,
 inu = void 0,
 j = void 0,
 k = void 0,
 nd = void 0,
 nlast = void 0,
 nn = void 0,
 nuf = void 0,
 nw = void 0,
 nz = void 0;

 bry = new Array(3);
 cssr = new Array(3);
 csrr = new Array(3);
 cyr = new Array(2);
 cyi = new Array(2);

 zeror = 0.0;
 zeroi = 0.0;
 coner = 1.0;

 cipr = [1, 0, -1, 0];
 cipi = [0, 1, 0, -1];
 hpi = 1.57079632679489662;
 aic = 1.265512123484645396;


 var goToLabel = 0;
 mainExecutionLoop: while (true) {
 switch (goToLabel) {
 case 0:
 nz = 0;
 nd = n;
 nlast = 0;
 // c-----------------------------------------------------------------------
 // c computed values with exponents between alim and elim in mag-
 // c nitude are scaled to keep intermediate arithmetic on scale,
 // c exp(alim)=exp(elim)*tol
 // c-----------------------------------------------------------------------
 cscl = 1.0 / tol;
 crsc = tol;
 cssr[0] = cscl;
 cssr[1] = coner;
 cssr[2] = crsc;
 csrr[0] = crsc;
 csrr[1] = coner;
 csrr[2] = cscl;
 bry[0] = 1.0e+3 * (0, _d1mach.d1mach)(1) / tol;
 // c-----------------------------------------------------------------------
 // c zn is in the right half plane after rotation by ci or -ci
 // c-----------------------------------------------------------------------
 znr = zi;
 zni = -zr;
 zbr = zr;
 zbi = zi;
 cidi = -coner;
 inu = Math.trunc(fnu);
 ang = hpi * (fnu - inu);
 c2r = Math.cos(ang);
 c2i = Math.sin(ang);
 car = c2r;
 sar = c2i;
 index = inu + n - 1;
 index = index % 4 + 1;
 str = c2r * cipr(index) - c2i * cipi(index);
 c2i = c2r * cipi(index) + c2i * cipr(index);
 c2r = str;
 if (zi > 0.0) {
 goToLabel = 10;break;
 }
 znr = -znr;
 zbi = -zbi;
 cidi = -cidi;
 c2i = -c2i;
 case 10:
 // c-----------------------------------------------------------------------
 // c check for underflow and overflow on first member
 // c-----------------------------------------------------------------------
 fn = Math.max(fnu, 1.0);

 var _zunhj = (0, _zunhj5.zunhj)(znr, zni, fn, 1, tol);

 var _zunhj2 = _slicedToArray(_zunhj, 12);

 phir = _zunhj2[0];
 phii = _zunhj2[1];
 argr = _zunhj2[2];
 argi = _zunhj2[3];
 zeta1r = _zunhj2[4];
 zeta1i = _zunhj2[5];
 zeta2r = _zunhj2[6];
 zeta2i = _zunhj2[7];
 asumr = _zunhj2[8];
 asumi = _zunhj2[9];
 bsumr = _zunhj2[10];
 bsumi = _zunhj2[11];

 if (kode === 1) {
 goToLabel = 20;break;
 }
 str = zbr + zeta2r;
 sti = zbi + zeta2i;
 rast = fn / (0, _zabs.azabs)(str, sti);
 str = str * rast * rast;
 sti = -sti * rast * rast;
 s1r = -zeta1r + str;
 s1i = -zeta1i + sti;
 goToLabel = 30;break;
 case 20:
 s1r = -zeta1r + zeta2r;
 s1i = -zeta1i + zeta2i;
 case 30:
 rs1 = s1r;
 if (Math.abs(rs1) > elim) {
 goToLabel = 150;break;
 }
 case 40:
 nn = Math.min(2, nd);
 // do 90 i=1,nn
 for (i = 1; i <= nn; i++) {
 fn = fnu + (nd - i);

 var _zunhj3 = (0, _zunhj5.zunhj)(znr, zni, fn, 0, tol);

 var _zunhj4 = _slicedToArray(_zunhj3, 12);

 phir = _zunhj4[0];
 phii = _zunhj4[1];
 argr = _zunhj4[2];
 argi = _zunhj4[3];
 zeta1r = _zunhj4[4];
 zeta1i = _zunhj4[5];
 zeta2r = _zunhj4[6];
 zeta2i = _zunhj4[7];
 asumr = _zunhj4[8];
 asumi = _zunhj4[9];
 bsumr = _zunhj4[10];
 bsumi = _zunhj4[11];

 if (kode === 1) {
 s1r = -zeta1r + zeta2r;
 s1i = -zeta1i + zeta2i;
 } else {
 str = zbr + zeta2r;
 sti = zbi + zeta2i;
 rast = fn / (0, _zabs.azabs)(str, sti);
 str = str * rast * rast;
 sti = -sti * rast * rast;
 s1r = -zeta1r + str;
 s1i = -zeta1i + sti + Math.abs(zi);
 }
 // c-----------------------------------------------------------------------
 // c test for underflow and overflow
 // c-----------------------------------------------------------------------
 rs1 = s1r;
 if (Math.abs(rs1) > elim) {
 goToLabel = 120;break;
 }
 if (i === 1) iflag = 2;
 if (Math.abs(rs1) < alim) {//
 // go to 70
 } else {
 // c-----------------------------------------------------------------------
 // c refine test and scale
 // c-----------------------------------------------------------------------
 // c-----------------------------------------------------------------------
 aphi = (0, _zabs.azabs)(phir, phii);
 aarg = (0, _zabs.azabs)(argr, argi);
 rs1 = rs1 + Math.log(aphi) - 0.25 * Math.log(aarg) - aic;
 if (Math.abs(rs1) > elim) {
 goToLabel = 120;break;
 }
 if (i === 1) iflag = 1;
 if (rs1 < 0.0) {
 // go to 70
 } else {
 if (i === 1) iflag = 3;
 }
 }
 // 70 continue
 // c-----------------------------------------------------------------------
 // c scale s1 to keep intermediate arithmetic on scale near
 // c exponent extremes
 // c-----------------------------------------------------------------------

 var _zairy = (0, _zairy5.zairy)(argr, argi, 0, 2);

 var _zairy2 = _slicedToArray(_zairy, 2);

 air = _zairy2[0];
 aii = _zairy2[1];

 var _zairy3 = (0, _zairy5.zairy)(argr, argi, 1, 2);

 var _zairy4 = _slicedToArray(_zairy3, 2);

 dair = _zairy4[0];
 daii = _zairy4[1];

 str = dair * bsumr - daii * bsumi;
 sti = dair * bsumi + daii * bsumr;
 str = str + (air * asumr - aii * asumi);
 sti = sti + (air * asumi + aii * asumr);
 s2r = phir * str - phii * sti;
 s2i = phir * sti + phii * str;
 str = Math.exp(s1r) * cssr[iflag - 1];
 s1r = str * Math.cos(s1i);
 s1i = str * Math.sin(s1i);
 str = s2r * s1r - s2i * s1i;
 s2i = s2r * s1i + s2i * s1r;
 s2r = str;
 if (iflag !== 1) {
 // go to 80
 } else {
 nw = (0, _zuchk.zuchk)(s2r, s2i, bry[0], tol);
 if (nw !== 0) {
 goToLabel = 120;break;
 }
 }
 // 80 continue
 if (zi <= 0.0) s2i = -s2i;
 str = s2r * c2r - s2i * c2i;
 s2i = s2r * c2i + s2i * c2r;
 s2r = str;
 cyr[i - 1] = s2r;
 cyi[i - 1] = s2i;
 j = nd - i + 1;
 yr[j - 1] = s2r * csrr[iflag - 1];
 yi[j - 1] = s2i * csrr[iflag - 1];
 str = -c2i * cidi;
 c2i = c2r * cidi;
 c2r = str;
 }
 if (goToLabel === 120) {
 break;
 }
 // 90 continue
 if (nd <= 2) {
 goToLabel = 110;break;
 }
 raz = 1.0 / (0, _zabs.azabs)(zr, zi);
 str = zr * raz;
 sti = -zi * raz;
 rzr = (str + str) * raz;
 rzi = (sti + sti) * raz;
 bry[1] = 1.0 / bry[0];
 bry[2] = (0, _d1mach.d1mach)(2);
 s1r = cyr[0];
 s1i = cyi[0];
 s2r = cyr[1];
 s2i = cyi[1];
 c1r = csrr[iflag - 1];
 ascle = bry[iflag - 1];
 k = nd - 2;
 fn = k;
 // do 100 i=3,nd
 for (i = 3; i <= nd; i++) {
 c2r = s2r;
 c2i = s2i;
 s2r = s1r + (fnu + fn) * (rzr * c2r - rzi * c2i);
 s2i = s1i + (fnu + fn) * (rzr * c2i + rzi * c2r);
 s1r = c2r;
 s1i = c2i;
 c2r = s2r * c1r;
 c2i = s2i * c1r;
 yr[k - 1] = c2r;
 yi[k - 1] = c2i;
 k = k - 1;
 fn = fn - 1.0;
 if (iflag >= 3) continue;
 str = Math.abs(c2r);
 sti = Math.abs(c2i);
 c2m = Math.max(str, sti);
 if (c2m <= ascle) continue;
 iflag = iflag + 1;
 ascle = bry[iflag - 1];
 s1r = s1r * c1r;
 s1i = s1i * c1r;
 s2r = c2r;
 s2i = c2i;
 s1r = s1r * cssr[iflag - 1];
 s1i = s1i * cssr[iflag - 1];
 s2r = s2r * cssr[iflag - 1];
 s2i = s2i * cssr[iflag - 1];
 c1r = csrr[iflag - 1];
 }
 // 100 continue
 case 110:
 break mainExecutionLoop;
 case 120:
 if (rs1 > 0.0) {
 goToLabel = 140;break;
 }
 // c-----------------------------------------------------------------------
 // c set underflow and update parameters
 // c-----------------------------------------------------------------------
 yr[nd - 1] = zeror;
 yi[nd - 1] = zeroi;
 nz = nz + 1;
 nd = nd - 1;
 if (nd === 0) {
 goToLabel = 110;break;
 }
 nuf = (0, _zuoik.zuoik)(zr, zi, fnu, kode, 1, nd, yr, yi, tol, elim, alim);
 if (nuf < 0) {
 goToLabel = 140;break;
 }
 nd = nd - nuf;
 nz = nz + nuf;
 if (nd === 0) {
 goToLabel = 110;break;
 }
 fn = fnu + (nd - 1);
 if (fn < fnul) {
 goToLabel = 130;break;
 }
 // c fn = cidi
 // c j = nuf + 1
 // c k = mod(j,4) + 1
 // c s1r = cipr(k)
 // c s1i = cipi(k)
 // c if (fn < 0.0) s1i = -s1i
 // c str = c2r*s1r - c2i*s1i
 // c c2i = c2r*s1i + c2i*s1r
 // c c2r = str
 index = inu + nd - 1;
 index = index % 4 + 1;
 c2r = car * cipr[index - 1] - sar * cipi[index - 1];
 c2i = car * cipi[index - 1] + sar * cipr[index - 1];
 if (zi <= 0.0) c2i = -c2i;
 goToLabel = 40;break;
 case 130:
 nlast = nd;
 break mainExecutionLoop;
 case 140:
 nz = -1;
 break mainExecutionLoop;
 case 150:
 if (rs1 > 0.0) {
 goToLabel = 140;break;
 }
 nz = n;
 // do 160 i=1,n
 for (i = 1; i <= n; i++) {
 yr[i - 1] = zeror;
 yi[i - 1] = zeroi;
 }
 // 160 continue
 default:
 break mainExecutionLoop;
 }
 }

 return [nz, nlast];
}
},{"../../utils/fortran-utils/d1mach.js":91,"./zabs.js":11,"./zairy.js":14,"./zuchk.js":37,"./zunhj.js":38,"./zuoik.js":44}],41:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});

var _slicedToArray = function () { function sliceIterator(arr, i) { var _arr = []; var _n = true; var _d = false; var _e = undefined; try { for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) { _arr.push(_s.value); if (i && _arr.length === i) break; } } catch (err) { _d = true; _e = err; } finally { try { if (!_n && _i["return"]) _i["return"](); } finally { if (_d) throw _e; } } return _arr; } return function (arr, i) { if (Array.isArray(arr)) { return arr; } else if (Symbol.iterator in Object(arr)) { return sliceIterator(arr, i); } else { throw new TypeError("Invalid attempt to destructure non-iterable instance"); } }; }(); // SUBROUTINE ZUNIK(ZRR, ZRI, FNU, IKFLG, IPMTR, TOL, INIT, PHIR,
// * PHII, ZETA1R, ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI)
// ***BEGIN PROLOGUE ZUNIK
// ***REFER TO ZBESI,ZBESK
//
// ZUNIK COMPUTES PARAMETERS FOR THE UNIFORM ASYMPTOTIC
// EXPANSIONS OF THE I AND K FUNCTIONS ON IKFLG= 1 OR 2
// RESPECTIVELY BY
//
// W(FNU,ZR) = PHI*EXP(ZETA)*SUM
//
// WHERE ZETA=-ZETA1 + ZETA2 OR
// ZETA1 - ZETA2
//
// THE FIRST CALL MUST HAVE INIT=0. SUBSEQUENT CALLS WITH THE
// SAME ZR AND FNU WILL RETURN THE I OR K FUNCTION ON IKFLG=
// 1 OR 2 WITH NO CHANGE IN INIT. CWRK IS A COMPLEX WORK
// ARRAY. IPMTR=0 COMPUTES ALL PARAMETERS. IPMTR=1 COMPUTES PHI,
// ZETA1,ZETA2.
//
// ***ROUTINES CALLED ZDIV,AZLOG,AZSQRT,D1MACH
// ***END PROLOGUE ZUNIK
// COMPLEX CFN,CON,CONE,CRFN,CWRK,CZERO,PHI,S,SR,SUM,T,T2,ZETA1,
// *ZETA2,ZN,ZR


exports.zunik = zunik;

var _zdiv7 = require('./zdiv.js');

var _zlog = require('./zlog.js');

var _zsqrt = require('./zsqrt.js');

var _d1mach = require('../../utils/fortran-utils/d1mach.js');

function zunik(zrr, zri, fnu, ikflg, ipmtr, tol, init) {
 var ac = void 0,
 c = void 0,
 con = void 0,
 conei = void 0,
 coner = void 0,
 crfni = void 0,
 crfnr = void 0,
 cwrki = void 0,
 cwrkr = void 0,
 phii = void 0,
 phir = void 0,
 rfn = void 0,
 si = void 0,
 sr = void 0,
 sri = void 0,
 srr = void 0,
 sti = void 0,
 str = void 0,
 sumi = void 0,
 sumr = void 0,
 test = void 0,
 ti = void 0,
 tr = void 0,
 t2i = void 0,
 t2r = void 0,
 zeroi = void 0,
 zeror = void 0,
 zeta1i = void 0,
 zeta1r = void 0,
 zeta2i = void 0,
 zeta2r = void 0,
 zni = void 0,
 znr = void 0,
 i = void 0,
 j = void 0,
 k = void 0,
 l = void 0;
 cwrkr = new Array(16);
 cwrki = new Array(16);
 zeror = 0;
 zeroi = 0;
 coner = 1;
 conei = 0;
 con = [3.98942280401432678e-01, 1.25331413731550025e+00];
 c = [1.00000000000000000e+00, -2.08333333333333333e-01, 1.25000000000000000e-01, 3.34201388888888889e-01, -4.01041666666666667e-01, 7.03125000000000000e-02, -1.02581259645061728e+00, 1.84646267361111111e+00, -8.91210937500000000e-01, 7.32421875000000000e-02, 4.66958442342624743e+00, -1.12070026162229938e+01, 8.78912353515625000e+00, -2.36408691406250000e+00, 1.12152099609375000e-01, -2.82120725582002449e+01, 8.46362176746007346e+01, -9.18182415432400174e+01, 4.25349987453884549e+01, -7.36879435947963170e+00, 2.27108001708984375e-01, 2.12570130039217123e+02, -7.65252468141181642e+02, 1.05999045252799988e+03, -6.99579627376132541e+02, 2.18190511744211590e+02, -2.64914304869515555e+01, 5.72501420974731445e-01, -1.91945766231840700e+03, 8.06172218173730938e+03, -1.35865500064341374e+04, 1.16553933368645332e+04, -5.30564697861340311e+03, 1.20090291321635246e+03, -1.08090919788394656e+02, 1.72772750258445740e+00, 2.02042913309661486e+04, -9.69805983886375135e+04, 1.92547001232531532e+05, -2.03400177280415534e+05, 1.22200464983017460e+05, -4.11926549688975513e+04, 7.10951430248936372e+03, -4.93915304773088012e+02, 6.07404200127348304e+00, -2.42919187900551333e+05, 1.31176361466297720e+06, -2.99801591853810675e+06, 3.76327129765640400e+06, -2.81356322658653411e+06, 1.26836527332162478e+06, -3.31645172484563578e+05, 4.52187689813627263e+04, -2.49983048181120962e+03, 2.43805296995560639e+01, 3.28446985307203782e+06, -1.97068191184322269e+07, 5.09526024926646422e+07, -7.41051482115326577e+07, 6.63445122747290267e+07, -3.75671766607633513e+07, 1.32887671664218183e+07, -2.78561812808645469e+06, 3.08186404612662398e+05, -1.38860897537170405e+04, 1.10017140269246738e+02, -4.93292536645099620e+07, 3.25573074185765749e+08, -9.39462359681578403e+08, 1.55359689957058006e+09, -1.62108055210833708e+09, 1.10684281682301447e+09, -4.95889784275030309e+08, 1.42062907797533095e+08, -2.44740627257387285e+07, 2.24376817792244943e+06, -8.40054336030240853e+04, 5.51335896122020586e+02, 8.14789096118312115e+08, -5.86648149205184723e+09, 1.86882075092958249e+10, -3.46320433881587779e+10, 4.12801855797539740e+10, -3.30265997498007231e+10, 1.79542137311556001e+10, -6.56329379261928433e+09, 1.55927986487925751e+09, -2.25105661889415278e+08, 1.73951075539781645e+07, -5.49842327572288687e+05, 3.03809051092238427e+03, -1.46792612476956167e+10, 1.14498237732025810e+11, -3.99096175224466498e+11, 8.19218669548577329e+11, -1.09837515608122331e+12, 1.00815810686538209e+12, -6.45364869245376503e+11, 2.87900649906150589e+11, -8.78670721780232657e+10, 1.76347306068349694e+10, -2.16716498322379509e+09, 1.43157876718888981e+08, -3.87183344257261262e+06, 1.82577554742931747e+04, 2.86464035717679043e+11, -2.40629790002850396e+12, 9.10934118523989896e+12, -2.05168994109344374e+13, 3.05651255199353206e+13, -3.16670885847851584e+13, 2.33483640445818409e+13, -1.23204913055982872e+13, 4.61272578084913197e+12, -1.19655288019618160e+12, 2.05914503232410016e+11, -2.18229277575292237e+10, 1.24700929351271032e+09, -2.91883881222208134e+07, 1.18838426256783253e+05];
 if (init !== 0) {
 // go to 40
 } else {
 // c-----------------------------------------------------------------------
 // c initialize all variables
 // c-----------------------------------------------------------------------
 rfn = 1.0 / fnu;
 // c-----------------------------------------------------------------------
 // c overflow test (zr/fnu too small)
 // c-----------------------------------------------------------------------
 test = (0, _d1mach.d1mach)(1) * 1.0e3;
 ac = fnu * test;
 if (Math.abs(zrr) > ac || Math.abs(zri) > ac) {
 // go to 15
 } else {
 zeta1r = 2.0 * Math.abs(Math.log(test)) + fnu;
 zeta1i = 0.0;
 zeta2r = fnu;
 zeta2i = 0.0;
 phir = 1.0;
 phii = 0.0;
 sumr = sumi = null;
 return [phir, phii, zeta1r, zeta1i, zeta2r, zeta2i, sumr, sumi, cwrkr, cwrki];
 }
 // 15 continue
 tr = zrr * rfn;
 ti = zri * rfn;
 sr = coner + (tr * tr - ti * ti);
 si = conei + (tr * ti + ti * tr);

 var _azsqrt = (0, _zsqrt.azsqrt)(sr, si);

 var _azsqrt2 = _slicedToArray(_azsqrt, 2);

 srr = _azsqrt2[0];
 sri = _azsqrt2[1];

 str = coner + srr;
 sti = conei + sri;

 var _zdiv = (0, _zdiv7.zdiv)(str, sti, tr, ti);

 var _zdiv2 = _slicedToArray(_zdiv, 2);

 znr = _zdiv2[0];
 zni = _zdiv2[1];

 var _azlog = (0, _zlog.azlog)(znr, zni);

 var _azlog2 = _slicedToArray(_azlog, 2);

 str = _azlog2[0];
 sti = _azlog2[1];

 zeta1r = fnu * str;
 zeta1i = fnu * sti;
 zeta2r = fnu * srr;
 zeta2i = fnu * sri;

 var _zdiv3 = (0, _zdiv7.zdiv)(coner, conei, srr, sri);

 var _zdiv4 = _slicedToArray(_zdiv3, 2);

 tr = _zdiv4[0];
 ti = _zdiv4[1];

 srr = tr * rfn;
 sri = ti * rfn;

 var _azsqrt3 = (0, _zsqrt.azsqrt)(srr, sri);

 var _azsqrt4 = _slicedToArray(_azsqrt3, 2);

 cwrkr[15] = _azsqrt4[0];
 cwrki[15] = _azsqrt4[1];

 phir = cwrkr[15] * con[ikflg - 1];
 phii = cwrki[15] * con[ikflg - 1];
 if (ipmtr !== 0) {
 return [phir, phii, zeta1r, zeta1i, zeta2r, zeta2i, sumr, sumi, cwrkr, cwrki];
 }

 var _zdiv5 = (0, _zdiv7.zdiv)(coner, conei, sr, si);

 var _zdiv6 = _slicedToArray(_zdiv5, 2);

 t2r = _zdiv6[0];
 t2i = _zdiv6[1];

 cwrkr[0] = coner;
 cwrki[0] = conei;
 crfnr = coner;
 crfni = conei;
 ac = 1.0;
 l = 1;
 // do 20 k=2,15
 for (k = 2; k <= 15; k++) {
 sr = zeror;
 si = zeroi;
 // do 10 j=1,k
 for (j = 1; j <= k; j++) {
 l = l + 1;
 str = sr * t2r - si * t2i + c[l - 1];
 si = sr * t2i + si * t2r;
 sr = str;
 }
 // 10 continue
 str = crfnr * srr - crfni * sri;
 crfni = crfnr * sri + crfni * srr;
 crfnr = str;
 cwrkr[k - 1] = crfnr * sr - crfni * si;
 cwrki[k - 1] = crfnr * si + crfni * sr;
 ac = ac * rfn;
 test = Math.abs(cwrkr[k - 1]) + Math.abs(cwrki[k - 1]);
 if (ac < tol && test < tol) {
 break;
 }
 }
 // 20 continue
 if (k === 16) k = 15; // if loop maxed out, set to 15, last index of loop; }
 // 30 continue
 init = k;
 }
 // 40 continue
 if (ikflg === 2) {
 // go to 60
 } else {
 // c-----------------------------------------------------------------------
 // c compute sum for the i function
 // c-----------------------------------------------------------------------
 sr = zeror;
 si = zeroi;
 // do 50 i=1,init
 for (i = 1; i <= init; i++) {
 sr = sr + cwrkr[i - 1];
 si = si + cwrki[i - 1];
 }
 // 50 continue
 sumr = sr;
 sumi = si;
 phir = cwrkr[15] * con[0];
 phii = cwrki[15] * con[0];
 return [phir, phii, zeta1r, zeta1i, zeta2r, zeta2i, sumr, sumi, cwrkr, cwrki];
 }
 // 60 continue
 // c-----------------------------------------------------------------------
 // c compute sum for the k function
 // c-----------------------------------------------------------------------
 sr = zeror;
 si = zeroi;
 tr = coner;
 // do 70 i=1,init
 for (i = 1; i <= init; i++) {
 sr = sr + tr * cwrkr(i - 1);
 si = si + tr * cwrki(i - 1);
 tr = -tr;
 }
 // 70 continue
 sumr = sr;
 sumi = si;
 phir = cwrkr[15] * con[1];
 phii = cwrki[15] * con[1];
 return [phir, phii, zeta1r, zeta1i, zeta2r, zeta2i, sumr, sumi, cwrkr, cwrki];
}
},{"../../utils/fortran-utils/d1mach.js":91,"./zdiv.js":26,"./zlog.js":29,"./zsqrt.js":36}],42:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});

var _slicedToArray = function () { function sliceIterator(arr, i) { var _arr = []; var _n = true; var _d = false; var _e = undefined; try { for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) { _arr.push(_s.value); if (i && _arr.length === i) break; } } catch (err) { _d = true; _e = err; } finally { try { if (!_n && _i["return"]) _i["return"](); } finally { if (_d) throw _e; } } return _arr; } return function (arr, i) { if (Array.isArray(arr)) { return arr; } else if (Symbol.iterator in Object(arr)) { return sliceIterator(arr, i); } else { throw new TypeError("Invalid attempt to destructure non-iterable instance"); } }; }(); /* eslint-disable no-labels */
/* eslint-disable no-fallthrough */
// SUBROUTINE ZUNK1(ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, TOL, ELIM,
// * ALIM)
// ***BEGIN PROLOGUE ZUNK1
// ***REFER TO ZBESK
//
// ZUNK1 COMPUTES K(FNU,Z) AND ITS ANALYTIC CONTINUATION FROM THE
// RIGHT HALF PLANE TO THE LEFT HALF PLANE BY MEANS OF THE
// UNIFORM ASYMPTOTIC EXPANSION.
// MR INDICATES THE DIRECTION OF ROTATION FOR ANALYTIC CONTINUATION.
// NZ=-1 MEANS AN OVERFLOW WILL OCCUR
//
// ***ROUTINES CALLED ZKSCL,ZS1S2,ZUCHK,ZUNIK,D1MACH,AZABS
// ***END PROLOGUE ZUNK1


exports.zunk1 = zunk1;

var _zs1s5 = require('./zs1s2.js');

var _zuchk = require('./zuchk.js');

var _zunik7 = require('./zunik.js');

var _zabs = require('./zabs.js');

var _d1mach = require('../../utils/fortran-utils/d1mach.js');

var _fortranHelpers = require('../../utils/fortranHelpers.js');

var ft = _interopRequireWildcard(_fortranHelpers);

function _interopRequireWildcard(obj) { if (obj && obj.__esModule) { return obj; } else { var newObj = {}; if (obj != null) { for (var key in obj) { if (Object.prototype.hasOwnProperty.call(obj, key)) newObj[key] = obj[key]; } } newObj.default = obj; return newObj; } }

function zunk1(zr, zi, fnu, kode, mr, n, yr, yi, tol, elim, alim) {
 var ang = void 0,
 aphi = void 0,
 asc = void 0,
 ascle = void 0,
 bry = void 0,
 cki = void 0,
 ckr = void 0,
 coner = void 0,
 crsc = void 0,
 cscl = void 0,
 csgni = void 0,
 cspni = void 0,
 cspnr = void 0,
 csr = void 0,
 csrr = void 0,
 cssr = void 0,
 cyi = void 0,
 cyr = void 0,
 c1i = void 0,
 c1r = void 0,
 c2i = void 0,
 c2m = void 0,
 c2r = void 0,
 fmr = void 0,
 fn = void 0,
 fnf = void 0,
 phidi = void 0,
 phidr = void 0,
 phii = void 0,
 phir = void 0,
 pi = void 0,
 rast = void 0,
 razr = void 0,
 rs1 = void 0,
 rzi = void 0,
 rzr = void 0,
 sgn = void 0,
 sti = void 0,
 str = void 0,
 sumdi = void 0,
 sumdr = void 0,
 sumi = void 0,
 sumr = void 0,
 s1i = void 0,
 s1r = void 0,
 s2i = void 0,
 s2r = void 0,
 zeroi = void 0,
 zeror = void 0,
 zeta1i = void 0,
 zeta1r = void 0,
 zeta2i = void 0,
 zeta2r = void 0,
 zet1di = void 0,
 zet1dr = void 0,
 zet2di = void 0,
 zet2dr = void 0,
 zri = void 0,
 zrr = void 0,
 i = void 0,
 ib = void 0,
 iflag = void 0,
 ifn = void 0,
 il = void 0,
 init = void 0,
 inu = void 0,
 iuf = void 0,
 k = void 0,
 kdflg = void 0,
 kflag = void 0,
 kk = void 0,
 nw = void 0,
 nz = void 0,
 initd = void 0,
 ic = void 0,
 ipard = void 0,
 j = void 0,
 f = void 0,
 m = void 0,
 zbr = void 0,
 zbi = void 0,
 yy = void 0;

 bry = new Array(3);
 init = new Array(2);
 sumr = new Array(2);
 sumi = new Array(2);
 zeta1r = new Array(2);
 zeta1i = new Array(2);
 zeta2r = new Array(2);
 zeta2i = new Array(2);
 cyr = new Array(2);
 cyi = new Array(2);
 cssr = new Array(3);
 csrr = new Array(3);
 phir = new Array(2);
 phii = new Array(2);

 // Init 2D arrs for cwrk
 var cwrkr = [];
 var cwrki = [];
 var iMax = 16;
 var jMax = 3;
 for (i = 0; i < iMax; i++) {
 cwrkr[i] = [];
 cwrki[i] = [];
 for (j = 0; j < jMax; j++) {
 f[i][j] = null;
 }
 }

 zeror = 0.0;
 zeroi = 0.0;
 coner = 1.0;

 pi = 3.14159265358979324;

 var goToLabel = 0;
 mainExecutionLoop: while (true) {
 switch (goToLabel) {
 case 0:
 kdflg = 1;
 nz = 0;
 // c-----------------------------------------------------------------------
 // c exp(-alim)=exp(-elim)/tol=approx. one precision greater than
 // c the underflow limit
 // c-----------------------------------------------------------------------
 cscl = 1.0 / tol;
 crsc = tol;
 cssr[0] = cscl;
 cssr[1] = coner;
 cssr[2] = crsc;
 csrr[0] = crsc;
 csrr[1] = coner;
 csrr[2] = cscl;
 bry[0] = 1.0e+3 * (0, _d1mach.d1mach)(1) / tol;
 bry[1] = 1.0 / bry[0];
 bry[2] = (0, _d1mach.d1mach)(2);
 zrr = zr;
 zri = zi;
 if (zr >= 0.0) {
 goToLabel = 10;break;
 }
 zrr = -zr;
 zri = -zi;
 case 10:
 j = 2;
 // do 70 i=1,n
 for (i = 1; i <= n; i++) {
 // c-----------------------------------------------------------------------
 // c j flip flops between 1 and 2 in j = 3 - j
 // c-----------------------------------------------------------------------
 j = 3 - j;
 fn = fnu + i - 1;
 init[j - 1] = 0;

 var _zunik = (0, _zunik7.zunik)(zrr, zri, fn, 2, 0, tol, init[j - 1]);

 var _zunik2 = _slicedToArray(_zunik, 10);

 phir[j - 1] = _zunik2[0];
 phii[j - 1] = _zunik2[1];
 zeta1r[j - 1] = _zunik2[2];
 zeta1i[j - 1] = _zunik2[3];
 zeta2r[j - 1] = _zunik2[4];
 zeta2i[j - 1] = _zunik2[5];
 sumr[j - 1] = _zunik2[6];
 sumi[j - 1] = _zunik2[7];
 cwrkr[0][j - 1] = _zunik2[8];
 cwrki[0][j - 1] = _zunik2[9];

 if (kode === 1) {
 s1r = zeta1r[j - 1] - zeta2r[j - 1];
 s1i = zeta1i[j - 1] - zeta2i[j - 1];
 } else {
 str = zrr + zeta2r[j - 1];
 sti = zri + zeta2i[j - 1];
 rast = fn / (0, _zabs.azabs)(str, sti);
 str = str * rast * rast;
 sti = -sti * rast * rast;
 s1r = zeta1r[j - 1] - str;
 s1i = zeta1i[j - 1] - sti;
 }
 rs1 = s1r;
 // c-----------------------------------------------------------------------
 // c test for underflow and overflow
 // c-----------------------------------------------------------------------
 if (Math.abs(rs1) > elim) goToLabel = 60;
 if (goToLabel !== 60) {
 if (kdflg === 1) kflag = 2;
 if (Math.abs(rs1) < alim) goToLabel = 40;
 if (goToLabel !== 40) {
 // c-----------------------------------------------------------------------
 // c refine test and scale
 // c-----------------------------------------------------------------------
 aphi = (0, _zabs.azabs)(phir[j - 1], phii[j - 1]);
 rs1 = rs1 + Math.log(aphi);
 if (Math.abs(rs1) > elim) goToLabel = 60;
 if (goToLabel !== 60) {
 if (kdflg === 1) kflag = 1;
 if (rs1 < 0.0) {
 // go to 50
 } else {
 if (kdflg === 1) kflag = 3;
 }
 }
 }
 // 40 continue
 if (goToLabel < 60) {
 // c-----------------------------------------------------------------------
 // c scale s1 to keep intermediate arithmetic on scale near
 // c exponent extremes
 // c-----------------------------------------------------------------------
 s2r = phir[j - 1] * sumr[j - 1] - phii[j - 1] * sumi[j - 1];
 s2i = phir[j - 1] * sumi[j - 1] + phii[j - 1] * sumr[j - 1];
 str = Math.exp(s1r) * cssr[kflag - 1];
 s1r = str * Math.cos(s1i);
 s1i = str * Math.sin(s1i);
 str = s2r * s1r - s2i * s1i;
 s2i = s1r * s2i + s2r * s1i;
 s2r = str;
 if (kflag !== 1) {
 // go to 50
 } else {
 nw = (0, _zuchk.zuchk)(s2r, s2i, bry[0], tol);
 if (nw !== 0) goToLabel = 60;
 }
 // 50 continue
 if (goToLabel < 60) {
 cyr[kdflg - 1] = s2r;
 cyi[kdflg - 1] = s2i;
 yr[i - 1] = s2r * csrr[kflag - 1];
 yi[i - 1] = s2i * csrr[kflag - 1];
 if (kdflg === 2) {
 goToLabel = 75;break;
 }
 kdflg = 2;
 break;
 }
 }
 }
 // 60 continue
 if (rs1 > 0.0) {
 goToLabel = 300;break;
 }
 // c-----------------------------------------------------------------------
 // c for zr < 0.0, the i function to be added will overflow
 // c-----------------------------------------------------------------------
 if (zr < 0.0) {
 goToLabel = 300;break;
 }
 kdflg = 1;
 yr[i - 1] = zeror;
 yi[i - 1] = zeroi;
 nz = nz + 1;
 if (i === 1) break;
 if (yr[i - 2] === zeror && yi[i - 2] === zeroi) break;
 yr[i - 2] = zeror;
 yi[i - 2] = zeroi;
 nz = nz + 1;
 }
 // 70 continue
 if (goToLabel < 70) {
 i = n;
 } else {
 break;
 }
 case 75:
 razr = 1.0 / (0, _zabs.azabs)(zrr, zri);
 str = zrr * razr;
 sti = -zri * razr;
 rzr = (str + str) * razr;
 rzi = (sti + sti) * razr;
 ckr = fn * rzr;
 cki = fn * rzi;
 ib = i + 1;
 if (n < ib) {
 goToLabel = 160;break;
 }
 // c-----------------------------------------------------------------------
 // c test last member for underflow and overflow. set sequence to zero
 // c on underflow.
 // c-----------------------------------------------------------------------
 fn = fnu + (n - 1);
 ipard = 1;
 if (mr !== 0) ipard = 0;
 initd = 0;

 var _zunik3 = (0, _zunik7.zunik)(zrr, zri, fn, 2, ipard, tol, initd);

 var _zunik4 = _slicedToArray(_zunik3, 10);

 phidr = _zunik4[0];
 phidi = _zunik4[1];
 zet1dr = _zunik4[2];
 zet1di = _zunik4[3];
 zet2dr = _zunik4[4];
 zet2di = _zunik4[5];
 sumdr = _zunik4[6];
 sumdi = _zunik4[7];
 cwrkr[0][2] = _zunik4[8];
 cwrki[0][2] = _zunik4[9];

 if (kode === 1) {
 goToLabel = 80;break;
 }
 str = zrr + zet2dr;
 sti = zri + zet2di;
 rast = fn / (0, _zabs.azabs)(str, sti);
 str = str * rast * rast;
 sti = -sti * rast * rast;
 s1r = zet1dr - str;
 s1i = zet1di - sti;
 goToLabel = 90;break;
 case 80:
 s1r = zet1dr - zet2dr;
 s1i = zet1di - zet2di;
 case 90:
 rs1 = s1r;
 if (Math.abs(rs1) > elim) {
 goToLabel = 95;break;
 }
 if (Math.abs(rs1) < alim) {
 goToLabel = 100;break;
 }
 // c----------------------------------------------------------------------------
 // c refine estimate and test
 // c-------------------------------------------------------------------------
 aphi = (0, _zabs.azabs)(phidr, phidi);
 rs1 = rs1 + Math.log(aphi);
 if (Math.abs(rs1) < elim) {
 goToLabel = 100;break;
 }
 case 95:
 if (Math.abs(rs1) > 0.0) {
 goToLabel = 300;break;
 }
 // c-----------------------------------------------------------------------
 // c for zr < 0.0, the i function to be added will overflow
 // c-----------------------------------------------------------------------
 if (zr < 0.0) {
 goToLabel = 300;break;
 }
 nz = n;
 // do 96 i=1,n
 for (i = 1; i <= n; i++) {
 yr[i - 1] = zeror;
 yi[i - 1] = zeroi;
 }
 // 96 continue
 break mainExecutionLoop;
 // c---------------------------------------------------------------------------
 // c forward recur for remainder of the sequence
 // c----------------------------------------------------------------------------
 case 100:
 s1r = cyr[0];
 s1i = cyi[0];
 s2r = cyr[1];
 s2i = cyi[1];
 c1r = csrr[kflag - 1];
 ascle = bry[kflag - 1];
 // do 120 i=ib,n
 for (i = ib; i <= n; i++) {
 c2r = s2r;
 c2i = s2i;
 s2r = ckr * c2r - cki * c2i + s1r;
 s2i = ckr * c2i + cki * c2r + s1i;
 s1r = c2r;
 s1i = c2i;
 ckr = ckr + rzr;
 cki = cki + rzi;
 c2r = s2r * c1r;
 c2i = s2i * c1r;
 yr[i - 1] = c2r;
 yi[i - 1] = c2i;
 if (kflag >= 3) break;
 str = Math.abs(c2r);
 sti = Math.abs(c2i);
 c2m = Math.max(str, sti);
 if (c2m <= ascle) break;
 kflag = kflag + 1;
 ascle = bry[kflag - 1];
 s1r = s1r * c1r;
 s1i = s1i * c1r;
 s2r = c2r;
 s2i = c2i;
 s1r = s1r * cssr[kflag - 1];
 s1i = s1i * cssr[kflag - 1];
 s2r = s2r * cssr[kflag - 1];
 s2i = s2i * cssr[kflag - 1];
 c1r = csrr[kflag - 1];
 }
 // 120 continue
 case 160:
 if (mr === 0) break mainExecutionLoop;
 // c-----------------------------------------------------------------------
 // c analytic continuation for re(z) < 0.0
 // c-----------------------------------------------------------------------
 nz = 0;
 fmr = mr;
 sgn = -ft.sign(pi, fmr);
 // c-----------------------------------------------------------------------
 // c cspn and csgn are coeff of k and i functions resp.
 // c-----------------------------------------------------------------------
 csgni = sgn;
 inu = Math.trunc(fnu);
 fnf = fnu - inu;
 ifn = inu + n - 1;
 ang = fnf * sgn;
 cspnr = Math.cos(ang);
 cspni = Math.sin(ang);
 if (ifn % 2 === 0) {
 goToLabel = 170;break;
 }
 cspnr = -cspnr;
 cspni = -cspni;
 case 170:
 asc = bry[0];
 iuf = 0;
 kk = n;
 kdflg = 1;
 ib = ib - 1;
 ic = ib - 1;
 // This loop 270 is such a hairy mess, I apologize if anyone has to come
 // back in here again. yeesh.
 // do 270 k=1,n
 loop270: for (k = 1; k <= n; k++) {
 fn = fnu + (kk - 1);
 // c-----------------------------------------------------------------------
 // c logic to sort out cases whose parameters were set for the k
 // c function above
 // c technically the below is correct but it's gross - KC
 // c-----------------------------------------------------------------------
 m = 3;
 if (n > 2 && kk === n && ib < n) {
 // do nothing
 } else if (n > 2 && (kk === ib || kk === ic) || n <= 2) {
 // 172 continue
 initd = init[j - 1];
 phidr = phir[j - 1];
 phidi = phii[j - 1];
 zet1dr = zeta1r[j - 1];
 zet1di = zeta1i[j - 1];
 zet2dr = zeta2r[j - 1];
 zet2di = zeta2i[j - 1];
 sumdr = sumr[j - 1];
 sumdi = sumi[j - 1];
 m = j;
 j = 3 - j;
 } else if (n > 2) {
 initd = 0;
 }
 // 180 continue

 var _zunik5 = (0, _zunik7.zunik)(zrr, zri, fn, 1, 0, tol, initd);

 var _zunik6 = _slicedToArray(_zunik5, 10);

 phidr = _zunik6[0];
 phidi = _zunik6[1];
 zet1dr = _zunik6[2];
 zet1di = _zunik6[3];
 zet2dr = _zunik6[4];
 zet2di = _zunik6[5];
 sumdr = _zunik6[6];
 sumdi = _zunik6[7];
 cwrkr[0][m - 1] = _zunik6[8];
 cwrki[0][m - 1] = _zunik6[9];

 if (kode === 1) {
 s1r = -zet1dr + zet2dr;
 s1i = -zet1di + zet2di;
 } else {
 str = zbr + zet2dr;
 sti = zbi + zet2di;
 rast = fn / (0, _zabs.azabs)(str, sti);
 str = str * rast * rast;
 sti = -sti * rast * rast;
 s1r = -zet1dr + str;
 s1i = -zet1di + sti;
 }
 // c-----------------------------------------------------------------------
 // c test for underflow and overflow
 // c-----------------------------------------------------------------------
 rs1 = s1r;
 if (Math.abs(rs1) > elim) {
 // go to 260
 if (rs1 > 0.0) {
 goToLabel = 300;break;
 }
 s2r = zeror;
 s2i = zeroi;
 // go to 230
 } else {
 if (kdflg === 1) iflag = 2;
 if (Math.abs(rs1) < alim) {
 // go to 220
 } else {
 // c-----------------------------------------------------------------------
 // c refine test and scale
 // c-----------------------------------------------------------------------
 aphi = (0, _zabs.azabs)(phidr, phidi);
 rs1 = rs1 + Math.log(aphi);
 if (Math.abs(rs1) > elim) {
 // go to 260
 if (rs1 > 0.0) {
 goToLabel = 300;break;
 }
 s2r = zeror;
 s2i = zeroi;
 goToLabel = 230;
 } else {
 if (kdflg === 1) iflag = 1;
 if (rs1 < 0.0) {
 // go to 220
 } else {
 if (kdflg === 1) iflag = 3;
 }
 }
 }
 // 220 continue
 if (goToLabel < 230) {
 str = phidr * sumdr - phidi * sumdi;
 sti = phidr * sumdi + phidi * sumdr;
 s2r = -csgni * sti;
 s2i = csgni * str;
 str = Math.exp(s1r) * cssr[iflag - 1];
 s1r = str * Math.cos(s1i);
 s1i = str * Math.sin(s1i);
 str = s2r * s1r - s2i * s1i;
 s2i = s2r * s1i + s2i * s1r;
 s2r = str;
 if (iflag !== 1) {
 // go to 230
 } else {
 nw = (0, _zuchk.zuchk)(s2r, s2i, bry[0], tol);
 if (nw === 0) {
 // go to 230
 } else {
 s2r = zeror;
 s2i = zeroi;
 }
 } // iflag !== 1
 } // if should run 220?
 } // if should skip straight to 230
 // 230 continue
 if (yy <= 0.0) s2i = -s2i;
 cyr[kdflg - 1] = s2r;
 cyi[kdflg - 1] = s2i;
 c2r = s2r;
 c2i = s2i;
 s2r = s2r * csrr[iflag - 1];
 s2i = s2i * csrr[iflag - 1];
 // c-----------------------------------------------------------------------
 // c add i and k functions, k sequence in y(i), i=1,n
 // c-----------------------------------------------------------------------
 s1r = yr[kk - 1];
 s1i = yi[kk - 1];
 if (kode === 1) {
 // go to 250
 } else {
 var _zs1s = (0, _zs1s5.zs1s2)(zrr, zri, s1r, s1i, s2r, s2i, asc, alim, iuf);

 var _zs1s2 = _slicedToArray(_zs1s, 5);

 s1r = _zs1s2[0];
 s1i = _zs1s2[1];
 s2r = _zs1s2[2];
 s2i = _zs1s2[3];
 nw = _zs1s2[4];

 nz = nz + nw;
 }
 // 250 continue
 yr[kk - 1] = s1r * cspnr - s1i * cspni + s2r;
 yi[kk - 1] = s1r * cspni + s1i * cspnr + s2i;
 kk = kk - 1;
 cspnr = -cspnr;
 cspni = -cspni;
 if (c2r !== 0.0 || c2i !== 0.0) {
 // go to 255
 } else {
 kdflg = 1;
 break;
 }
 // 255 continue
 if (kdflg === 2) {
 goToLabel = 295;break;
 }
 kdflg = 2;
 break;
 }
 // 270 continue
 k = n;
 case 275:
 il = n - k;
 if (il === 0) break mainExecutionLoop;
 // c-----------------------------------------------------------------------
 // c recur backward for remainder of i sequence and add in the
 // c k functions, scaling the i sequence during recurrence to keep
 // c intermediate arithmetic on scale near exponent extremes.
 // c-----------------------------------------------------------------------
 s1r = cyr[0];
 s1i = cyi[0];
 s2r = cyr[1];
 s2i = cyi[1];
 csr = csrr[iflag - 1];
 ascle = bry[iflag - 1];
 fn = inu + il;
 // do 290 i=1,il
 for (i = 1; i <= il; i++) {
 c2r = s2r;
 c2i = s2i;
 s2r = s1r + (fn + fnf) * (rzr * c2r - rzi * c2i);
 s2i = s1i + (fn + fnf) * (rzr * c2i + rzi * c2r);
 s1r = c2r;
 s1i = c2i;
 fn = fn - 1.0;
 c2r = s2r * csr;
 c2i = s2i * csr;
 ckr = c2r;
 cki = c2i;
 c1r = yr[kk - 1];
 c1i = yi[kk - 1];
 if (kode === 1) {
 // go to 280
 } else {
 var _zs1s3 = (0, _zs1s5.zs1s2)(zrr, zri, c1r, c1i, c2r, c2i, asc, alim, iuf);

 var _zs1s4 = _slicedToArray(_zs1s3, 5);

 c1r = _zs1s4[0];
 c1i = _zs1s4[1];
 c2r = _zs1s4[2];
 c2i = _zs1s4[3];
 nw = _zs1s4[4];

 nz = nz + nw;
 }
 // 280 continue
 yr[kk - 1] = c1r * cspnr - c1i * cspni + c2r;
 yi[kk - 1] = c1r * cspni + c1i * cspnr + c2i;
 kk = kk - 1;
 cspnr = -cspnr;
 cspni = -cspni;
 if (iflag >= 3) break;
 c2r = Math.abs(ckr);
 c2i = Math.abs(cki);
 c2m = Math.max(c2r, c2i);
 if (c2m <= ascle) break;
 iflag = iflag + 1;
 ascle = bry[iflag - 1];
 s1r = s1r * csr;
 s1i = s1i * csr;
 s2r = ckr;
 s2i = cki;
 s1r = s1r * cssr[iflag - 1];
 s1i = s1i * cssr[iflag - 1];
 s2r = s2r * cssr[iflag - 1];
 s2i = s2i * cssr[iflag - 1];
 csr = csrr[iflag - 1];
 }
 // 290 continue
 break mainExecutionLoop;
 case 300:
 nz = -1;
 default:
 break mainExecutionLoop;
 }
 }

 return nz;
}
},{"../../utils/fortran-utils/d1mach.js":91,"../../utils/fortranHelpers.js":93,"./zabs.js":11,"./zs1s2.js":33,"./zuchk.js":37,"./zunik.js":41}],43:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});

var _slicedToArray = function () { function sliceIterator(arr, i) { var _arr = []; var _n = true; var _d = false; var _e = undefined; try { for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) { _arr.push(_s.value); if (i && _arr.length === i) break; } } catch (err) { _d = true; _e = err; } finally { try { if (!_n && _i["return"]) _i["return"](); } finally { if (_d) throw _e; } } return _arr; } return function (arr, i) { if (Array.isArray(arr)) { return arr; } else if (Symbol.iterator in Object(arr)) { return sliceIterator(arr, i); } else { throw new TypeError("Invalid attempt to destructure non-iterable instance"); } }; }(); /* eslint-disable no-labels */
/* eslint-disable no-fallthrough */
// SUBROUTINE ZUNK2(ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, TOL, ELIM,
// * ALIM)
// ***BEGIN PROLOGUE ZUNK2
// ***REFER TO ZBESK
//
// ZUNK2 COMPUTES K(FNU,Z) AND ITS ANALYTIC CONTINUATION FROM THE
// RIGHT HALF PLANE TO THE LEFT HALF PLANE BY MEANS OF THE
// UNIFORM ASYMPTOTIC EXPANSIONS FOR H(KIND,FNU,ZN) AND J(FNU,ZN)
// WHERE ZN IS IN THE RIGHT HALF PLANE, KIND=(3-MR)/2, MR=+1 OR
// -1. HERE ZN=ZR*I OR -ZR*I WHERE ZR=Z IF Z IS IN THE RIGHT
// HALF PLANE OR ZR=-Z IF Z IS IN THE LEFT HALF PLANE. MR INDIC-
// ATES THE DIRECTION OF ROTATION FOR ANALYTIC CONTINUATION.
// NZ=-1 MEANS AN OVERFLOW WILL OCCUR
//
// ***ROUTINES CALLED ZAIRY,ZKSCL,ZS1S2,ZUCHK,ZUNHJ,D1MACH,AZABS
// ***END PROLOGUE ZUNK2


exports.zunk2 = zunk2;

var _zairy9 = require('./zairy.js');

var _zs1s5 = require('./zs1s2.js');

var _zuchk = require('./zuchk.js');

var _zunhj7 = require('./zunhj.js');

var _d1mach = require('../../utils/fortran-utils/d1mach.js');

var _zabs = require('./zabs.js');

var _fortranHelpers = require('../../utils/fortranHelpers.js');

var ft = _interopRequireWildcard(_fortranHelpers);

function _interopRequireWildcard(obj) { if (obj && obj.__esModule) { return obj; } else { var newObj = {}; if (obj != null) { for (var key in obj) { if (Object.prototype.hasOwnProperty.call(obj, key)) newObj[key] = obj[key]; } } newObj.default = obj; return newObj; } }

function zunk2(zr, zi, fnu, kode, mr, n, yr, yi, tol, elim, alim) {
 var aarg = void 0,
 aic = void 0,
 aii = void 0,
 air = void 0,
 ang = void 0,
 aphi = void 0,
 argdi = void 0,
 argdr = void 0,
 argi = void 0,
 argr = void 0,
 asc = void 0,
 ascle = void 0,
 asumdi = void 0,
 asumdr = void 0,
 asumi = void 0,
 asumr = void 0,
 bry = void 0,
 bsumdi = void 0,
 bsumdr = void 0,
 bsumi = void 0,
 bsumr = void 0,
 car = void 0,
 cipi = void 0,
 cipr = void 0,
 cki = void 0,
 ckr = void 0,
 coner = void 0,
 crsc = void 0,
 cr1i = void 0,
 cr1r = void 0,
 cr2i = void 0,
 cr2r = void 0,
 cscl = void 0,
 csgni = void 0,
 csi = void 0,
 cspni = void 0,
 cspnr = void 0,
 csr = void 0,
 csrr = void 0,
 cssr = void 0,
 cyi = void 0,
 cyr = void 0,
 c1i = void 0,
 c1r = void 0,
 c2i = void 0,
 c2m = void 0,
 c2r = void 0,
 daii = void 0,
 dair = void 0,
 fmr = void 0,
 fn = void 0,
 fnf = void 0,
 hpi = void 0,
 phidi = void 0,
 phidr = void 0,
 phii = void 0,
 phir = void 0,
 pi = void 0,
 pti = void 0,
 ptr = void 0,
 rast = void 0,
 razr = void 0,
 rs1 = void 0,
 rzi = void 0,
 rzr = void 0,
 sar = void 0,
 sgn = void 0,
 sti = void 0,
 str = void 0,
 s1i = void 0,
 s1r = void 0,
 s2i = void 0,
 s2r = void 0,
 yy = void 0,
 zbi = void 0,
 zbr = void 0,
 zeroi = void 0,
 zeror = void 0,
 zeta1i = void 0,
 zeta1r = void 0,
 zeta2i = void 0,
 zeta2r = void 0,
 zet1di = void 0,
 zet1dr = void 0,
 zet2di = void 0,
 zet2dr = void 0,
 zni = void 0,
 znr = void 0,
 zri = void 0,
 zrr = void 0,
 i = void 0,
 ib = void 0,
 iflag = void 0,
 ifn = void 0,
 il = void 0,
 index = void 0,
 inu = void 0,
 iuf = void 0,
 k = void 0,
 kdflg = void 0,
 kflag = void 0,
 kk = void 0,
 nw = void 0,
 nz = void 0,
 j = void 0,
 ipard = void 0,
 ic = void 0;
 bry = new Array(3);
 asumr = new Array(2);
 asumi = new Array(2);
 bsumr = new Array(2);
 bsumi = new Array(2);
 phir = new Array(2);
 phii = new Array(2);
 argr = new Array(2);
 argi = new Array(2);
 zeta1r = new Array(2);
 zeta1i = new Array(2);
 zeta2r = new Array(2);
 zeta2i = new Array(2);
 cyr = new Array(2);
 cyi = new Array(2);
 cssr = new Array(3);
 csrr = new Array(3);

 zeror = 0.0;
 zeroi = 0.0;
 coner = 1.0;
 cr1r = 1.0;
 cr1i = 1.73205080756887729;
 cr2r = -0.5;
 cr2i = -8.66025403784438647e-01;
 hpi = 1.57079632679489662;
 pi = 3.14159265358979324;
 aic = 1.26551212348464539;

 cipr = [1, 0, -1, 0];
 cipi = [0, -1, 0, 1];

 var goToLabel = 0;
 mainExecutionLoop: while (true) {
 switch (goToLabel) {
 case 0:
 kdflg = 1;
 nz = 0;
 // c-----------------------------------------------------------------------
 // c exp(-alim)=exp(-elim)/tol=approx. one precision greater than
 // c the underflow limit
 // c-----------------------------------------------------------------------
 cscl = 1.0 / tol;
 crsc = tol;
 cssr[0] = cscl;
 cssr[1] = coner;
 cssr[2] = crsc;
 csrr[0] = crsc;
 csrr[1] = coner;
 csrr[2] = cscl;
 bry[0] = 1.0e+3 * (0, _d1mach.d1mach)(1) / tol;
 bry[1] = 1.0 / bry[0];
 bry[2] = (0, _d1mach.d1mach)(2);
 zrr = zr;
 zri = zi;
 if (zr >= 0.0) {
 goToLabel = 10;break;
 }
 zrr = -zr;
 zri = -zi;
 case 10:
 yy = zri;
 znr = zri;
 zni = -zrr;
 zbr = zrr;
 zbi = zri;
 inu = Math.trunc(fnu);
 fnf = fnu - inu;
 ang = -hpi * fnf;
 car = Math.cos(ang);
 sar = Math.sin(ang);
 c2r = hpi * sar;
 c2i = -hpi * car;
 kk = inu % 4 + 1;
 str = c2r * cipr[kk - 1] - c2i * cipi[kk - 1];
 sti = c2r * cipi[kk - 1] + c2i * cipr[kk - 1];
 csr = cr1r * str - cr1i * sti;
 csi = cr1r * sti + cr1i * str;
 if (yy > 0.0) {
 goToLabel = 20;break;
 }
 znr = -znr;
 zbi = -zbi;
 case 20:
 // c-----------------------------------------------------------------------
 // c k(fnu,z) is computed from h(2,fnu,-i*z) where z is in the first
 // c quadrant. fourth quadrant values (yy <= 0.0) are computed by
 // c conjugation since the k function is real on the positive real axis
 // c-----------------------------------------------------------------------
 j = 2;
 // do 80 i=1,n
 for (i = 1; i <= n; i++) {
 // c-----------------------------------------------------------------------
 // c j flip flops between 1 and 2 in j = 3 - j
 // c-----------------------------------------------------------------------
 j = 3 - j;
 fn = fnu + i - 1;

 var _zunhj = (0, _zunhj7.zunhj)(znr, zni, fn, 0, tol);

 var _zunhj2 = _slicedToArray(_zunhj, 12);

 phir[j - 1] = _zunhj2[0];
 phii[j - 1] = _zunhj2[1];
 argr[j - 1] = _zunhj2[2];
 argi[j - 1] = _zunhj2[3];
 zeta1r[j - 1] = _zunhj2[4];
 zeta1i[j - 1] = _zunhj2[5];
 zeta2r[j - 1] = _zunhj2[6];
 zeta2i[j - 1] = _zunhj2[7];
 asumr[j - 1] = _zunhj2[8];
 asumi[j - 1] = _zunhj2[9];
 bsumr[j - 1] = _zunhj2[10];
 bsumi[j - 1] = _zunhj2[11];

 if (kode === 1) {
 s1r = zeta1r[j - 1] - zeta2r[j - 1];
 s1i = zeta1i[j - 1] - zeta2i[j - 1];
 } else {
 str = zbr + zeta2r[j - 1];
 sti = zbi + zeta2i[j - 1];
 rast = fn / (0, _zabs.azabs)(str, sti);
 str = str * rast * rast;
 sti = -sti * rast * rast;
 s1r = zeta1r[j - 1] - str;
 s1i = zeta1i[j - 1] - sti;
 }
 // c-----------------------------------------------------------------------
 // c test for underflow and overflow
 // c-----------------------------------------------------------------------
 rs1 = s1r;
 if (Math.abs(rs1) > elim) goToLabel = 70;
 if (goToLabel !== 70) {
 if (kdflg === 1) kflag = 2;
 if (Math.abs(rs1) < alim) goToLabel = 50;
 if (goToLabel !== 50) {
 // c-----------------------------------------------------------------------
 // c refine test and scale
 // c-----------------------------------------------------------------------
 aphi = (0, _zabs.azabs)(phir[j - 1], phii[j - 1]);
 aarg = (0, _zabs.azabs)(argr[j - 1], argi[j - 1]);
 rs1 = rs1 + Math.log(aphi) - 0.25 * Math.log(aarg) - aic;
 if (Math.abs(rs1) > elim) goToLabel = 70;
 if (goToLabel !== 70) {
 if (kdflg === 1) kflag = 1;
 if (rs1 < 0.0) {
 // go to 50
 } else {
 if (kdflg === 1) kflag = 3;
 }
 }
 }
 // 50 continue
 if (goToLabel < 70) {
 // c-----------------------------------------------------------------------
 // c scale s1 to keep intermediate arithmetic on scale near
 // c exponent extremes
 // c-----------------------------------------------------------------------
 c2r = argr[j - 1] * cr2r - argi[j - 1] * cr2i;
 c2i = argr[j - 1] * cr2i + argi[j - 1] * cr2r;

 var _zairy = (0, _zairy9.zairy)(c2r, c2i, 0, 2);

 var _zairy2 = _slicedToArray(_zairy, 2);

 air = _zairy2[0];
 aii = _zairy2[1];

 var _zairy3 = (0, _zairy9.zairy)(c2r, c2i, 1, 2);

 var _zairy4 = _slicedToArray(_zairy3, 2);

 dair = _zairy4[0];
 daii = _zairy4[1];

 str = dair * bsumr[j - 1] - daii * bsumi[j - 1];
 sti = dair * bsumi[j - 1] + daii * bsumr[j - 1];
 ptr = str * cr2r - sti * cr2i;
 pti = str * cr2i + sti * cr2r;
 str = ptr + (air * asumr[j - 1] - aii * asumi[j - 1]);
 sti = pti + (air * asumi[j - 1] + aii * asumr[j - 1]);
 ptr = str * phir[j - 1] - sti * phii[j - 1];
 pti = str * phii[j - 1] + sti * phir[j - 1];
 s2r = ptr * csr - pti * csi;
 s2i = ptr * csi + pti * csr;
 str = Math.exp(s1r) * cssr[kflag - 1];
 s1r = str * Math.cos(s1i);
 s1i = str * Math.sin(s1i);
 str = s2r * s1r - s2i * s1i;
 s2i = s1r * s2i + s2r * s1i;
 s2r = str;
 if (kflag !== 1) {
 // go to 60
 } else {
 nw = (0, _zuchk.zuchk)(s2r, s2i, bry[0], tol);
 if (nw !== 0) goToLabel = 70;
 }
 // 60 continue
 if (goToLabel < 70) {
 if (yy <= 0.0) s2i = -s2i;
 cyr[kdflg - 1] = s2r;
 cyi[kdflg - 1] = s2i;
 yr[i - 1] = s2r * csrr[kflag - 1];
 yi[i - 1] = s2i * csrr[kflag - 1];
 str = csi;
 csi = -csr;
 csr = str;
 if (kdflg === 2) {
 goToLabel = 85;break;
 }
 kdflg = 2;
 break;
 }
 }
 }
 // 70 continue
 if (rs1 > 0.0) {
 goToLabel = 320;break;
 }
 // c-----------------------------------------------------------------------
 // c for zr < 0.0, the i function to be added will overflow
 // c-----------------------------------------------------------------------
 if (zr < 0.0) {
 goToLabel = 320;break;
 }
 kdflg = 1;
 yr[i - 1] = zeror;
 yi[i - 1] = zeroi;
 nz = nz + 1;
 str = csi;
 csi = -csr;
 csr = str;
 if (i === 1) break;
 if (yr[i - 2] === zeror && yi[i - 2] === zeroi) break;
 yr[i - 2] = zeror;
 yi[i - 2] = zeroi;
 nz = nz + 1;
 }
 // 80 continue
 if (goToLabel < 80) {
 i = n;
 } else {
 break;
 }
 case 85:
 razr = 1.0 / (0, _zabs.azabs)(zrr, zri);
 str = zrr * razr;
 sti = -zri * razr;
 rzr = (str + str) * razr;
 rzi = (sti + sti) * razr;
 ckr = fn * rzr;
 cki = fn * rzi;
 ib = i + 1;
 if (n < ib) {
 goToLabel = 180;break;
 }
 // c-----------------------------------------------------------------------
 // c test last member for underflow and overflow. set sequence to zero
 // c on underflow.
 // c-----------------------------------------------------------------------
 fn = fnu + (n - 1);
 ipard = 1;
 if (mr !== 0) {
 ipard = 0;

 var _zunhj3 = (0, _zunhj7.zunhj)(znr, zni, fn, ipard, tol);

 var _zunhj4 = _slicedToArray(_zunhj3, 12);

 phidr = _zunhj4[0];
 phidi = _zunhj4[1];
 argdr = _zunhj4[2];
 argdi = _zunhj4[3];
 zet1dr = _zunhj4[4];
 zet1di = _zunhj4[5];
 zet2dr = _zunhj4[6];
 zet2di = _zunhj4[7];
 asumdr = _zunhj4[8];
 asumdi = _zunhj4[9];
 bsumdr = _zunhj4[10];
 bsumdi = _zunhj4[11];
 }
 if (kode === 1) {
 goToLabel = 90;break;
 }
 str = zbr + zet2dr;
 sti = zbi + zet2di;
 rast = fn / (0, _zabs.azabs)(str, sti);
 str = str * rast * rast;
 sti = -sti * rast * rast;
 s1r = zet1dr - str;
 s1i = zet1di - sti;
 goToLabel = 100;break;
 case 90:
 s1r = zet1dr - zet2dr;
 s1i = zet1di - zet2di;
 case 100:
 rs1 = s1r;
 if (Math.abs(rs1) > elim) {
 goToLabel = 105;break;
 }
 if (Math.abs(rs1) < alim) {
 goToLabel = 120;break;
 }
 // c----------------------------------------------------------------------------
 // c refine estimate and test
 // c-------------------------------------------------------------------------
 aphi = (0, _zabs.azabs)(phidr, phidi);
 rs1 = rs1 + Math.log(aphi);
 if (Math.abs(rs1) < elim) {
 goToLabel = 120;break;
 }
 case 105:
 if (rs1 > 0.0) {
 goToLabel = 320;break;
 }
 // c-----------------------------------------------------------------------
 // c for zr < 0.0, the i function to be added will overflow
 // c-----------------------------------------------------------------------
 if (zr < 0.0) {
 goToLabel = 320;break;
 }
 nz = n;
 // do 106 i=1,n
 for (i = 1; i <= n; i++) {
 yr[i - 1] = zeror;
 yi[i - 1] = zeroi;
 }
 // 106 continue
 break mainExecutionLoop;
 case 120:
 s1r = cyr[0];
 s1i = cyi[0];
 s2r = cyr[1];
 s2i = cyi[1];
 c1r = csrr[kflag - 1];
 ascle = bry[kflag - 1];
 // do 130 i=ib,n
 for (i = ib; i <= n; i++) {
 c2r = s2r;
 c2i = s2i;
 s2r = ckr * c2r - cki * c2i + s1r;
 s2i = ckr * c2i + cki * c2r + s1i;
 s1r = c2r;
 s1i = c2i;
 ckr = ckr + rzr;
 cki = cki + rzi;
 c2r = s2r * c1r;
 c2i = s2i * c1r;
 yr[i - 1] = c2r;
 yi[i - 1] = c2i;
 if (kflag >= 3) break;
 str = Math.abs(c2r);
 sti = Math.abs(c2i);
 c2m = Math.max(str, sti);
 if (c2m <= ascle) break;
 kflag = kflag + 1;
 ascle = bry[kflag - 1];
 s1r = s1r * c1r;
 s1i = s1i * c1r;
 s2r = c2r;
 s2i = c2i;
 s1r = s1r * cssr[kflag - 1];
 s1i = s1i * cssr[kflag - 1];
 s2r = s2r * cssr[kflag - 1];
 s2i = s2i * cssr[kflag - 1];
 c1r = csrr[kflag - 1];
 }
 // 130 continue
 case 180:
 if (mr === 0) break mainExecutionLoop;
 // c-----------------------------------------------------------------------
 // c analytic continuation for re(z) < 0.0
 // c-----------------------------------------------------------------------
 nz = 0;
 fmr = mr;
 sgn = -ft.sign(pi, fmr);
 // c-----------------------------------------------------------------------
 // c cspn and csgn are coeff of k and i functions resp.
 // c-----------------------------------------------------------------------
 csgni = sgn;
 if (yy <= 0.0) csgni = -csgni;
 ifn = inu + n - 1;
 ang = fnf * sgn;
 cspnr = Math.cos(ang);
 cspni = Math.sin(ang);
 if (ifn % 2 === 0) {
 goToLabel = 190;break;
 }
 cspnr = -cspnr;
 cspni = -cspni;
 case 190:
 // c-----------------------------------------------------------------------
 // c cs=coeff of the j function to get the i function. i(fnu,z) is
 // c computed from exp(i*fnu*hpi)*j(fnu,-i*z) where z is in the first
 // c quadrant. fourth quadrant values (yy <= 0.0e0) are computed by
 // c conjugation since the i function is real on the positive real axis
 // c-----------------------------------------------------------------------
 csr = sar * csgni;
 csi = car * csgni;
 index = ifn % 4 + 1;
 c2r = cipr[index - 1];
 c2i = cipi[index - 1];
 str = csr * c2r + csi * c2i;
 csi = -csr * c2i + csi * c2r;
 csr = str;
 asc = bry[0];
 iuf = 0;
 kk = n;
 kdflg = 1;
 ib = ib - 1;
 ic = ib - 1;
 // This loop 290 is such a hairy mess, I apologize if anyone has to come
 // back in here again. yeesh.
 // do 290 k=1,n
 for (k = 1; k <= n; k++) {
 fn = fnu + (kk - 1);
 // c-----------------------------------------------------------------------
 // c logic to sort out cases whose parameters were set for the k
 // c function above
 // c technically the below is correct but it's gross - KC
 // c-----------------------------------------------------------------------
 if (n > 2 && kk === n && ib < n) {
 // do nothing
 } else if (n > 2 && (kk === ib || kk === ic) || n <= 2) {
 // 172 continue
 phidr = phir[j - 1];
 phidi = phii[j - 1];
 argdr = argr[j - 1];
 argdi = argi[j - 1];
 zet1dr = zeta1r[j - 1];
 zet1di = zeta1i[j - 1];
 zet2dr = zeta2r[j - 1];
 zet2di = zeta2i[j - 1];
 asumdr = asumr[j - 1];
 asumdi = asumi[j - 1];
 bsumdr = bsumr[j - 1];
 bsumdi = bsumi[j - 1];
 j = 3 - j;
 } else if (n > 2) {
 var _zunhj5 = (0, _zunhj7.zunhj)(znr, zni, fn, 0, tol);

 var _zunhj6 = _slicedToArray(_zunhj5, 12);

 phidr = _zunhj6[0];
 phidi = _zunhj6[1];
 argdr = _zunhj6[2];
 argdi = _zunhj6[3];
 zet1dr = _zunhj6[4];
 zet1di = _zunhj6[5];
 zet2dr = _zunhj6[6];
 zet2di = _zunhj6[7];
 asumdr = _zunhj6[8];
 asumdi = _zunhj6[9];
 bsumdr = _zunhj6[10];
 bsumdi = _zunhj6[11];
 }
 // 210 continue
 if (kode === 1) {
 s1r = -zet1dr + zet2dr;
 s1i = -zet1di + zet2di;
 } else {
 str = zbr + zet2dr;
 sti = zbi + zet2di;
 rast = fn / (0, _zabs.azabs)(str, sti);
 str = str * rast * rast;
 sti = -sti * rast * rast;
 s1r = -zet1dr + str;
 s1i = -zet1di + sti;
 }
 // c-----------------------------------------------------------------------
 // c test for underflow and overflow
 // c-----------------------------------------------------------------------
 rs1 = s1r;
 if (Math.abs(rs1) > elim) {
 // go to 280
 if (rs1 > 0.0) {
 goToLabel = 320;break;
 }
 s2r = zeror;
 s2i = zeroi;
 // go to 250
 } else {
 if (kdflg === 1) iflag = 2;
 if (Math.abs(rs1) < alim) {
 // go to 240
 } else {
 // c-----------------------------------------------------------------------
 // c refine test and scale
 // c-----------------------------------------------------------------------
 aphi = (0, _zabs.azabs)(phidr, phidi);
 aarg = (0, _zabs.azabs)(argdr, argdi);
 rs1 = rs1 + Math.log(aphi) - 0.25 * Math.log(aarg) - aic;
 if (Math.abs(rs1) > elim) {
 // go to 280
 if (rs1 > 0.0) {
 goToLabel = 320;break;
 }
 s2r = zeror;
 s2i = zeroi;
 goToLabel = 250;
 } else {
 if (kdflg === 1) iflag = 1;
 if (rs1 < 0.0) {
 // go to 240
 } else {
 if (kdflg === 1) iflag = 3;
 }
 }
 }
 // 240 continue
 if (goToLabel < 250) {
 var _zairy5 = (0, _zairy9.zairy)(argdr, argdi, 0, 2);

 var _zairy6 = _slicedToArray(_zairy5, 2);

 air = _zairy6[0];
 aii = _zairy6[1];

 var _zairy7 = (0, _zairy9.zairy)(argdr, argdi, 1, 2);

 var _zairy8 = _slicedToArray(_zairy7, 2);

 dair = _zairy8[0];
 daii = _zairy8[1];

 str = dair * bsumdr - daii * bsumdi;
 sti = dair * bsumdi + daii * bsumdr;
 str = str + (air * asumdr - aii * asumdi);
 sti = sti + (air * asumdi + aii * asumdr);
 ptr = str * phidr - sti * phidi;
 pti = str * phidi + sti * phidr;
 s2r = ptr * csr - pti * csi;
 s2i = ptr * csi + pti * csr;
 str = Math.exp(s1r) * cssr[iflag - 1];
 s1r = str * Math.cos(s1i);
 s1i = str * Math.sin(s1i);
 str = s2r * s1r - s2i * s1i;
 s2i = s2r * s1i + s2i * s1r;
 s2r = str;
 if (iflag !== 1) {
 // go to 250
 } else {
 nw = (0, _zuchk.zuchk)(s2r, s2i, bry[0], tol);
 if (nw === 0) {
 // go to 250
 } else {
 s2r = zeror;
 s2i = zeroi;
 }
 } // iflag !== 1
 } // if should run 240?
 } // if should skip straight to 250
 // 250 continue
 if (yy <= 0.0) s2i = -s2i;
 cyr[kdflg - 1] = s2r;
 cyi[kdflg - 1] = s2i;
 c2r = s2r;
 c2i = s2i;
 s2r = s2r * csrr[iflag - 1];
 s2i = s2i * csrr[iflag - 1];
 // c-----------------------------------------------------------------------
 // c add i and k functions, k sequence in y(i), i=1,n
 // c-----------------------------------------------------------------------
 s1r = yr[kk - 1];
 s1i = yi[kk - 1];
 if (kode === 1) {
 // go to 270
 } else {
 var _zs1s = (0, _zs1s5.zs1s2)(zrr, zri, s1r, s1i, s2r, s2i, asc, alim, iuf);

 var _zs1s2 = _slicedToArray(_zs1s, 5);

 s1r = _zs1s2[0];
 s1i = _zs1s2[1];
 s2r = _zs1s2[2];
 s2i = _zs1s2[3];
 nw = _zs1s2[4];

 nz = nz + nw;
 }
 // 270 continue
 yr[kk - 1] = s1r * cspnr - s1i * cspni + s2r;
 yi[kk - 1] = s1r * cspni + s1i * cspnr + s2i;
 kk = kk - 1;
 cspnr = -cspnr;
 cspni = -cspni;
 str = csi;
 csi = -csr;
 csr = str;
 if (c2r !== 0.0 || c2i !== 0.0) {
 // go to 275
 } else {
 kdflg = 1;
 break;
 }
 // 275 continue
 if (kdflg === 2) {
 goToLabel = 295;break;
 }
 kdflg = 2;
 break;
 }
 // 290 continue
 k = n;
 case 295:
 il = n - k;
 if (il === 0) break mainExecutionLoop;
 // c-----------------------------------------------------------------------
 // c recur backward for remainder of i sequence and add in the
 // c k functions, scaling the i sequence during recurrence to keep
 // c intermediate arithmetic on scale near exponent extremes.
 // c-----------------------------------------------------------------------
 s1r = cyr[0];
 s1i = cyi[0];
 s2r = cyr[1];
 s2i = cyi[1];
 csr = csrr[iflag - 1];
 ascle = bry[iflag - 1];
 fn = inu + il;
 // do 310 i=1,il
 for (i = 1; i <= il; i++) {
 c2r = s2r;
 c2i = s2i;
 s2r = s1r + (fn + fnf) * (rzr * c2r - rzi * c2i);
 s2i = s1i + (fn + fnf) * (rzr * c2i + rzi * c2r);
 s1r = c2r;
 s1i = c2i;
 fn = fn - 1.0;
 c2r = s2r * csr;
 c2i = s2i * csr;
 ckr = c2r;
 cki = c2i;
 c1r = yr[kk - 1];
 c1i = yi[kk - 1];
 if (kode === 1) {
 // go to 300
 } else {
 var _zs1s3 = (0, _zs1s5.zs1s2)(zrr, zri, c1r, c1i, c2r, c2i, asc, alim, iuf);

 var _zs1s4 = _slicedToArray(_zs1s3, 5);

 c1r = _zs1s4[0];
 c1i = _zs1s4[1];
 c2r = _zs1s4[2];
 c2i = _zs1s4[3];
 nw = _zs1s4[4];

 nz = nz + nw;
 }
 // 300 continue
 yr[kk - 1] = c1r * cspnr - c1i * cspni + c2r;
 yi[kk - 1] = c1r * cspni + c1i * cspnr + c2i;
 kk = kk - 1;
 cspnr = -cspnr;
 cspni = -cspni;
 if (iflag >= 3) break;
 c2r = Math.abs(ckr);
 c2i = Math.abs(cki);
 c2m = Math.max(c2r, c2i);
 if (c2m <= ascle) break;
 iflag = iflag + 1;
 ascle = bry[iflag - 1];
 s1r = s1r * csr;
 s1i = s1i * csr;
 s2r = ckr;
 s2i = cki;
 s1r = s1r * cssr[iflag - 1];
 s1i = s1i * cssr[iflag - 1];
 s2r = s2r * cssr[iflag - 1];
 s2i = s2i * cssr[iflag - 1];
 csr = csrr[iflag - 1];
 }
 // 310 continue
 break mainExecutionLoop;
 case 320:
 nz = -1;
 default:
 break mainExecutionLoop;
 }
 }

 return nz;
}
},{"../../utils/fortran-utils/d1mach.js":91,"../../utils/fortranHelpers.js":93,"./zabs.js":11,"./zairy.js":14,"./zs1s2.js":33,"./zuchk.js":37,"./zunhj.js":38}],44:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});

var _slicedToArray = function () { function sliceIterator(arr, i) { var _arr = []; var _n = true; var _d = false; var _e = undefined; try { for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) { _arr.push(_s.value); if (i && _arr.length === i) break; } } catch (err) { _d = true; _e = err; } finally { try { if (!_n && _i["return"]) _i["return"](); } finally { if (_d) throw _e; } } return _arr; } return function (arr, i) { if (Array.isArray(arr)) { return arr; } else if (Symbol.iterator in Object(arr)) { return sliceIterator(arr, i); } else { throw new TypeError("Invalid attempt to destructure non-iterable instance"); } }; }(); /* eslint-disable no-labels */
/* eslint-disable no-fallthrough */
// SUBROUTINE ZUOIK(ZR, ZI, FNU, KODE, IKFLG, N, YR, YI, NUF, TOL,
// * ELIM, ALIM)
// ***BEGIN PROLOGUE ZUOIK
// ***REFER TO ZBESI,ZBESK,ZBESH
//
// ZUOIK COMPUTES THE LEADING TERMS OF THE UNIFORM ASYMPTOTIC
// EXPANSIONS FOR THE I AND K FUNCTIONS AND COMPARES THEM
// (IN LOGARITHMIC FORM) TO ALIM AND ELIM FOR OVER AND UNDERFLOW
// WHERE ALIM.LT.ELIM. IF THE MAGNITUDE, BASED ON THE LEADING
// EXPONENTIAL, IS LESS THAN ALIM OR GREATER THAN -ALIM, THEN
// THE RESULT IS ON SCALE. IF NOT, THEN A REFINED TEST USING OTHER
// MULTIPLIERS (IN LOGARITHMIC FORM) IS MADE BASED ON ELIM. HERE
// EXP(-ELIM)=SMALLEST MACHINE NUMBER*1.0E+3 AND EXP(-ALIM)=
// EXP(-ELIM)/TOL
//
// IKFLG=1 MEANS THE I SEQUENCE IS TESTED
// =2 MEANS THE K SEQUENCE IS TESTED
// NUF = 0 MEANS THE LAST MEMBER OF THE SEQUENCE IS ON SCALE
// =-1 MEANS AN OVERFLOW WOULD OCCUR
// IKFLG=1 AND NUF.GT.0 MEANS THE LAST NUF Y VALUES WERE SET TO ZERO
// THE FIRST N-NUF VALUES MUST BE SET BY ANOTHER ROUTINE
// IKFLG=2 AND NUF.EQ.N MEANS ALL Y VALUES WERE SET TO ZERO
// IKFLG=2 AND 0.LT.NUF.LT.N NOT CONSIDERED. Y MUST BE SET BY
// ANOTHER ROUTINE
//
// ***ROUTINES CALLED ZUCHK,ZUNHJ,ZUNIK,D1MACH,AZABS,AZLOG
// ***END PROLOGUE ZUOIK
// COMPLEX ARG,ASUM,BSUM,CWRK,CZ,CZERO,PHI,SUM,Y,Z,ZB,ZETA1,ZETA2,ZN,
// *ZR


exports.zuoik = zuoik;

var _d1mach = require('../../utils/fortran-utils/d1mach.js');

var _zabs = require('./zabs.js');

var _zlog = require('./zlog.js');

var _zuchk = require('./zuchk.js');

var _zunhj5 = require('./zunhj.js');

var _zunik5 = require('./zunik.js');

function zuoik(zr, zi, fnu, kode, ikflg, n, yr, yi, tol, elim, alim) {
 var aarg = void 0,
 aic = void 0,
 aphi = void 0,
 argi = void 0,
 argr = void 0,
 ascle = void 0,
 ax = void 0,
 ay = void 0,
 czi = void 0,
 czr = void 0,
 fnn = void 0,
 gnn = void 0,
 gnu = void 0,
 phii = void 0,
 phir = void 0,
 rcz = void 0,
 str = void 0,
 sti = void 0,
 zbi = void 0,
 zbr = void 0,
 zeroi = void 0,
 zeror = void 0,
 zeta1i = void 0,
 zeta1r = void 0,
 zeta2i = void 0,
 zeta2r = void 0,
 zni = void 0,
 znr = void 0,
 zri = void 0,
 zrr = void 0,
 i = void 0,
 iform = void 0,
 init = void 0,
 nn = void 0,
 nuf = void 0,
 nw = void 0;

 zeror = 0;
 zeroi = 0;

 aic = 1.265512123484645396;

 var goToLabel = 0;
 mainExecutionLoop: while (true) {
 switch (goToLabel) {
 case 0:
 nuf = 0;
 nn = n;
 zrr = zr;
 zri = zi;
 if (zr >= 0.0) {
 goToLabel = 10;break;
 }
 zrr = -zr;
 zri = -zi;
 case 10:
 zbr = zrr;
 zbi = zri;
 ax = Math.abs(zr) * 1.7321;
 ay = Math.abs(zi);
 iform = 1;
 if (ay > ax) iform = 2;
 gnu = Math.max(fnu, 1.0);
 if (ikflg === 1) {
 goToLabel = 20;break;
 }
 fnn = nn;
 gnn = fnu + fnn - 1.0;
 gnu = Math.max(gnn, fnn);
 case 20:
 // c-----------------------------------------------------------------------
 // c only the magnitude of arg and phi are needed along with the
 // c real parts of zeta1, zeta2 and zb. no attempt is made to get
 // c the sign of the imaginary part correct.
 // c-----------------------------------------------------------------------
 if (iform === 2) {
 goToLabel = 30;break;
 }
 init = 0;

 var _zunik = (0, _zunik5.zunik)(zrr, zri, gnu, ikflg, 1, tol, init);

 var _zunik2 = _slicedToArray(_zunik, 6);

 phir = _zunik2[0];
 phii = _zunik2[1];
 zeta1r = _zunik2[2];
 zeta1i = _zunik2[3];
 zeta2r = _zunik2[4];
 zeta2i = _zunik2[5];

 czr = -zeta1r + zeta2r;
 czi = -zeta1i + zeta2i;
 goToLabel = 50;break;
 case 30:
 znr = zri;
 zni = -zrr;
 if (zi > 0.0) {
 goToLabel = 40;break;
 }
 znr = -znr;
 case 40:
 var _zunhj = (0, _zunhj5.zunhj)(znr, zni, gnu, 1, tol);

 var _zunhj2 = _slicedToArray(_zunhj, 8);

 phir = _zunhj2[0];
 phii = _zunhj2[1];
 argr = _zunhj2[2];
 argi = _zunhj2[3];
 zeta1r = _zunhj2[4];
 zeta1i = _zunhj2[5];
 zeta2r = _zunhj2[6];
 zeta2i = _zunhj2[7];

 czr = -zeta1r + zeta2r;
 czi = -zeta1i + zeta2i;
 aarg = (0, _zabs.azabs)(argr, argi);
 case 50:
 if (kode === 1) {
 goToLabel = 60;break;
 }
 czr = czr - zbr;
 czi = czi - zbi;
 case 60:
 if (ikflg === 1) {
 goToLabel = 70;break;
 }
 czr = -czr;
 czi = -czi;
 case 70:
 aphi = (0, _zabs.azabs)(phir, phii);
 rcz = czr;
 // c-----------------------------------------------------------------------
 // c overflow test
 // c-----------------------------------------------------------------------
 if (rcz > elim) {
 goToLabel = 210;break;
 }
 if (rcz < alim) {
 goToLabel = 80;break;
 }
 rcz = rcz + Math.log(aphi);
 if (iform === 2) rcz = rcz - 0.25 * Math.log(aarg) - aic;
 if (rcz > elim) {
 goToLabel = 210;break;
 }
 goToLabel = 130;break;
 case 80:
 // c-----------------------------------------------------------------------
 // c underflow test
 // c-----------------------------------------------------------------------
 if (rcz < -elim) {
 goToLabel = 90;break;
 }
 if (rcz > -alim) {
 goToLabel = 130;break;
 }
 rcz = rcz + Math.log(aphi);
 if (iform === 2) rcz = rcz - 0.25 * Math.log(aarg) - aic;
 if (rcz > -elim) {
 goToLabel = 110;break;
 }
 case 90:
 // do 100 i=1,nn
 for (i = 1; i <= nn; i++) {
 yr[i - 1] = zeror;
 yi[i - 1] = zeroi;
 }
 // 100 dcontinue
 nuf = nn;
 break mainExecutionLoop;
 case 110:
 ascle = 1.0e+3 * (0, _d1mach.d1mach)(1) / tol;

 var _azlog = (0, _zlog.azlog)(phir, phii);

 var _azlog2 = _slicedToArray(_azlog, 2);

 str = _azlog2[0];
 sti = _azlog2[1];

 czr = czr + str;
 czi = czi + sti;
 if (iform === 1) {
 goToLabel = 120;break;
 }

 var _azlog3 = (0, _zlog.azlog)(argr, argi);

 var _azlog4 = _slicedToArray(_azlog3, 2);

 str = _azlog4[0];
 sti = _azlog4[1];

 czr = czr - 0.25 * str - aic;
 czi = czi - 0.25 * sti;
 case 120:
 ax = Math.exp(rcz) / tol;
 ay = czi;
 czr = ax * Math.cos(ay);
 czi = ax * Math.sin(ay);
 nw = (0, _zuchk.zuchk)(czr, czi, ascle, tol);
 if (nw !== 0) {
 goToLabel = 90;break;
 }
 case 130:
 if (ikflg === 2) break mainExecutionLoop;
 if (n === 1) break mainExecutionLoop;
 // c-----------------------------------------------------------------------
 // c set underflows on i sequence
 // c-----------------------------------------------------------------------
 case 140:
 gnu = fnu + (nn - 1);
 if (iform === 2) {
 goToLabel = 150;break;
 }
 init = 0;

 var _zunik3 = (0, _zunik5.zunik)(zrr, zri, gnu, ikflg, 1, tol, init);

 var _zunik4 = _slicedToArray(_zunik3, 6);

 phir = _zunik4[0];
 phii = _zunik4[1];
 zeta1r = _zunik4[2];
 zeta1i = _zunik4[3];
 zeta2r = _zunik4[4];
 zeta2i = _zunik4[5];

 czr = -zeta1r + zeta2r;
 czi = -zeta1i + zeta2i;
 goToLabel = 160;break;
 case 150:
 var _zunhj3 = (0, _zunhj5.zunhj)(znr, zni, gnu, 1, tol);

 var _zunhj4 = _slicedToArray(_zunhj3, 8);

 phir = _zunhj4[0];
 phii = _zunhj4[1];
 argr = _zunhj4[2];
 argi = _zunhj4[3];
 zeta1r = _zunhj4[4];
 zeta1i = _zunhj4[5];
 zeta2r = _zunhj4[6];
 zeta2i = _zunhj4[7];

 czr = -zeta1r + zeta2r;
 czi = -zeta1i + zeta2i;
 aarg = (0, _zabs.azabs)(argr, argi);
 case 160:
 if (kode === 1) {
 goToLabel = 170;break;
 }
 czr = czr - zbr;
 czi = czi - zbi;
 case 170:
 aphi = (0, _zabs.azabs)(phir, phii);
 rcz = czr;
 if (rcz < -elim) {
 goToLabel = 180;break;
 }
 if (rcz > -alim) break mainExecutionLoop;
 rcz = rcz + Math.log(aphi);
 if (iform === 2) rcz = rcz - 0.25 * Math.log(aarg) - aic;
 if (rcz > -elim) {
 goToLabel = 190;break;
 }
 case 180:
 yr[nn - 1] = zeror;
 yi[nn - 1] = zeroi;
 nn = nn - 1;
 nuf = nuf + 1;
 if (nn === 0) break mainExecutionLoop;
 goToLabel = 140;break;
 case 190:
 ascle = 1.0e+3 * (0, _d1mach.d1mach)(1) / tol;

 var _azlog5 = (0, _zlog.azlog)(phir, phii);

 var _azlog6 = _slicedToArray(_azlog5, 2);

 str = _azlog6[0];
 sti = _azlog6[1];

 czr = czr + str;
 czi = czi + sti;
 if (iform === 1) {
 goToLabel = 200;break;
 }

 var _azlog7 = (0, _zlog.azlog)(argr, argi);

 var _azlog8 = _slicedToArray(_azlog7, 2);

 str = _azlog8[0];
 sti = _azlog8[1];

 czr = czr - 0.25 * str - aic;
 czi = czi - 0.25 * sti;
 case 200:
 ax = Math.exp(rcz) / tol;
 ay = czi;
 czr = ax * Math.cos(ay);
 czi = ax * Math.sin(ay);
 nw = (0, _zuchk.zuchk)(czr, czi, ascle, tol);
 if (nw !== 0) {
 goToLabel = 180;break;
 }
 break mainExecutionLoop;
 case 210:
 nuf = -1;
 default:
 break mainExecutionLoop;
 }
 }

 return nuf;
}
},{"../../utils/fortran-utils/d1mach.js":91,"./zabs.js":11,"./zlog.js":29,"./zuchk.js":37,"./zunhj.js":38,"./zunik.js":41}],45:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.zwrsk = zwrsk;

var _d1mach = require('../../utils/fortran-utils/d1mach.js');

var _zbknu = require('./zbknu.js');

var _zrati = require('./zrati.js');

var _zabs = require('./zabs.js');

/* eslint-disable no-labels */
/* eslint-disable no-fallthrough */
// SUBROUTINE ZWRSK(ZRR, ZRI, FNU, KODE, N, YR, YI, NZ, CWR, CWI,
// * TOL, ELIM, ALIM)
// ***BEGIN PROLOGUE ZWRSK
// ***REFER TO ZBESI,ZBESK
//
// ZWRSK COMPUTES THE I BESSEL FUNCTION FOR RE(Z).GE.0.0 BY
// NORMALIZING THE I FUNCTION RATIOS FROM ZRATI BY THE WRONSKIAN
//
// ***ROUTINES CALLED D1MACH,ZBKNU,ZRATI,AZABS
// ***END PROLOGUE ZWRSK
function zwrsk(zrr, zri, fnu, kode, n, yr, yi, cwr, cwi, tol, elim, alim) {
 var act = void 0,
 acw = void 0,
 ascle = void 0,
 cinui = void 0,
 cinur = void 0,
 csclr = void 0,
 cti = void 0,
 ctr = void 0,
 c1i = void 0,
 c1r = void 0,
 c2i = void 0,
 c2r = void 0,
 pti = void 0,
 ptr = void 0,
 ract = void 0,
 sti = void 0,
 str = void 0,
 i = void 0,
 nw = void 0,
 nz = void 0;
 cwr = new Array(2);
 cwi = new Array(2);
 var goToLabel = 0;
 mainExecutionLoop: while (true) {
 switch (goToLabel) {
 case 0:
 // c-----------------------------------------------------------------------
 // c i(fnu+i-1,z) by backward recurrence for ratios
 // c y(i)=i(fnu+i,z)/i(fnu+i-1,z) from crati normalized by the
 // c wronskian with k(fnu,z) and k(fnu+1,z) from cbknu.
 // c-----------------------------------------------------------------------
 nz = 0;
 nw = (0, _zbknu.zbknu)(zrr, zri, fnu, kode, 2, cwr, cwi, tol, elim, alim);
 if (nw !== 0) {
 goToLabel = 50;break;
 }
 (0, _zrati.zrati)(zrr, zri, fnu, n, yr, yi, tol);
 // c-----------------------------------------------------------------------
 // c recur forward on i(fnu+1,z) = r(fnu,z)*i(fnu,z),
 // c r(fnu+j-1,z)=y(j), j=1,...,n
 // c-----------------------------------------------------------------------
 cinur = 1.0;
 cinui = 0.0;
 if (kode === 1) {
 goToLabel = 10;break;
 }
 cinur = Math.cos(zri);
 cinui = Math.sin(zri);
 case 10:
 // c-----------------------------------------------------------------------
 // c on low exponent machines the k functions can be close to both
 // c the under and overflow limits and the normalization must be
 // c scaled to prevent over or underflow. cuoik has determined that
 // c the result is on scale.
 // c-----------------------------------------------------------------------
 acw = (0, _zabs.azabs)(cwr[1], cwi[1]);
 ascle = 1.0e+3 * (0, _d1mach.d1mach)(1) / tol;
 csclr = 1.0;
 if (acw > ascle) {
 goToLabel = 20;break;
 }
 csclr = 1.0 / tol;
 goToLabel = 30;break;
 case 20:
 ascle = 1.0 / ascle;
 if (acw < ascle) {
 goToLabel = 30;break;
 }
 csclr = tol;
 case 30:
 c1r = cwr[0] * csclr;
 c1i = cwi[0] * csclr;
 c2r = cwr[1] * csclr;
 c2i = cwi[1] * csclr;
 str = yr[0];
 sti = yi[0];
 // c-----------------------------------------------------------------------
 // c cinu=cinu*(conjg(ct)/cabs(ct))*(1.0/cabs(ct) prevents
 // c under- or overflow prematurely by squaring cabs(ct)
 // c-----------------------------------------------------------------------
 ptr = str * c1r - sti * c1i;
 pti = str * c1i + sti * c1r;
 ptr = ptr + c2r;
 pti = pti + c2i;
 ctr = zrr * ptr - zri * pti;
 cti = zrr * pti + zri * ptr;
 act = (0, _zabs.azabs)(ctr, cti);
 ract = 1.0 / act;
 ctr = ctr * ract;
 cti = -cti * ract;
 ptr = cinur * ract;
 pti = cinui * ract;
 cinur = ptr * ctr - pti * cti;
 cinui = ptr * cti + pti * ctr;
 yr[0] = cinur * csclr;
 yi[0] = cinui * csclr;
 if (n === 1) break mainExecutionLoop;
 // do 40 i=2,n
 for (i = 2; i <= n; i++) {
 ptr = str * cinur - sti * cinui;
 cinui = str * cinui + sti * cinur;
 cinur = ptr;
 str = yr[i - 1];
 sti = yi[i - 1];
 yr[i - 1] = cinur * csclr;
 yi[i - 1] = cinui * csclr;
 }
 // 40 continue
 break mainExecutionLoop;
 case 50:
 nz = -1;
 if (nw === -2) nz = -2;
 default:
 break mainExecutionLoop;
 }
 }

 return nz;
}
},{"../../utils/fortran-utils/d1mach.js":91,"./zabs.js":11,"./zbknu.js":23,"./zrati.js":32}],46:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
function gammasgn(x) {
 var fx = void 0;

 if (isNaN(x)) {
 return x;
 }
 if (x > 0) {
 return 1.0;
 } else {
 fx = Math.floor(x);
 if (x - fx === 0.0) {
 return 0.0;
 } else if (fx % 2) {
 return -1.0;
 } else {
 return 1.0;
 }
 }
}

exports.gammasgn = gammasgn;
},{}],47:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.struveL = exports.struveH = undefined;

var _slicedToArray = function () { function sliceIterator(arr, i) { var _arr = []; var _n = true; var _d = false; var _e = undefined; try { for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) { _arr.push(_s.value); if (i && _arr.length === i) break; } } catch (err) { _d = true; _e = err; } finally { try { if (!_n && _i["return"]) _i["return"](); } finally { if (_d) throw _e; } } return _arr; } return function (arr, i) { if (Array.isArray(arr)) { return arr; } else if (Symbol.iterator in Object(arr)) { return sliceIterator(arr, i); } else { throw new TypeError("Invalid attempt to destructure non-iterable instance"); } }; }(); /*
 * Compute the Struve function.
 *
 * Notes
 * -----
 *
 * We use three expansions for the Struve function discussed in [1]:
 *
 * - power series
 * - expansion in Bessel functions
 * - asymptotic large-z expansion
 *
 * Rounding errors are estimated based on the largest terms in the sums.
 *
 * ``struve_convergence.js`` *will* plot the convergence regions of the different
 * expansions (someday).
 *
 * (i)
 *
 * Looking at the error in the asymptotic expansion, one finds that
 * it's not worth trying if z ~> 0.7 * v + 12 for v > 0.
 *
 * (ii)
 *
 * The Bessel function expansion tends to fail for |z| >~ |v| and is not tried
 * there.
 *
 * For Struve H it covers the quadrant v > z where the power series may fail to
 * produce reasonable results.
 *
 * (iii)
 *
 * The three expansions together cover for Struve H the region z > 0, v real.
 *
 * They also cover Struve L, except that some loss of precision may occur around
 * the transition region z ~ 0.7 |v|, v < 0, |v| >> 1 where the function changes
 * rapidly.
 *
 * (iv)
 *
 * The power series is evaluated in double-double precision. This fixes accuracy
 * issues in Struve H for |v| << |z| before the asymptotic expansion kicks in.
 * Moreover, it improves the Struve L behavior for negative v.
 *
 *
 * References
 * ----------
 * [1] NIST Digital Library of Mathematical Functions
 * https://dlmf.nist.gov/11
 */

/*
 * Copyright (C) 2013 Pauli Virtanen
 * Ported to ECMAScript 2018 by KC Erb
 * Copyright (C) 2018, Kings Distributed Systems
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * a. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 * b. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 * c. Neither the name of Enthought nor the names of the SciPy Developers
 * may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
 * OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 */


var _gammasgn = require('./gammasgn.js');

var _gamma = require('../cephes/gamma.js');

var _lgam = require('../cephes/gamma/lgam.js');

var _cephes = require('../cephes.js');

var _cephes2 = _interopRequireDefault(_cephes);

var _dd = require('../cephes/dd.js');

var _dd2 = _interopRequireDefault(_dd);

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

var MAXITER = 10000;
var SUM_EPS = 1e-16; /* be sure we are in the tail of the sum */
var SUM_TINY = 1e-100;
var GOOD_EPS = 1e-12;
var ACCEPTABLE_EPS = 1e-7;
var ACCEPTABLE_ATOL = 1e-300;

function struveH(v, z) {
 return struveHL(v, z, 1);
}

function struveL(v, z) {
 return struveHL(v, z, 0);
}

function struveHL(v, z, isH) {
 var value = void 0,
 err = void 0,
 tmp = void 0,
 n = void 0;
 value = new Float64Array(4);
 err = new Float64Array(4);

 if (z < 0) {
 n = Math.trunc(v);
 if (v === n) {
 tmp = n % 2 === 0 ? -1 : 1;
 return tmp * struveHL(v, -z, isH);
 } else {
 return NaN;
 }
 } else if (z === 0) {
 if (v < -1) {
 return (0, _gammasgn.gammasgn)(v + 1.5) * Infinity;
 } else if (v === -1) {
 return 2 / Math.sqrt(Math.PI) / (0, _gamma.gamma)(0.5);
 } else {
 return 0;
 }
 }

 n = -v - 0.5;
 if (n === -v - 0.5 && n > 0) {
 if (isH) {
 return (n % 2 === 0 ? 1 : -1) * _cephes2.default.jv(n + 0.5, z);
 } else {
 return _cephes2.default.iv(n + 0.5, z);
 }
 }

 /* Try the asymptotic expansion */
 if (z >= 0.7 * v + 12) {
 var _struveAsympLargeZ = struveAsympLargeZ(v, z, isH);

 var _struveAsympLargeZ2 = _slicedToArray(_struveAsympLargeZ, 2);

 value[0] = _struveAsympLargeZ2[0];
 err[0] = _struveAsympLargeZ2[1];

 if (err[0] < GOOD_EPS * Math.abs(value[0])) {
 return value[0];
 }
 } else {
 err[0] = Infinity;
 }

 /* Try power series */

 var _struvePowerSeries = struvePowerSeries(v, z, isH);

 var _struvePowerSeries2 = _slicedToArray(_struvePowerSeries, 2);

 value[1] = _struvePowerSeries2[0];
 err[1] = _struvePowerSeries2[1];

 if (err[1] < GOOD_EPS * Math.abs(value[1])) {
 return value[1];
 }

 /* Try bessel series */
 if (Math.abs(z) < Math.abs(v) + 20) {
 var _struveBesselSeries = struveBesselSeries(v, z, isH);

 var _struveBesselSeries2 = _slicedToArray(_struveBesselSeries, 2);

 value[2] = _struveBesselSeries2[0];
 err[2] = _struveBesselSeries2[1];

 if (err[2] < GOOD_EPS * Math.abs(value[2])) {
 return value[2];
 }
 } else {
 err[2] = Infinity;
 }

 /* Return the best of the three, if it is acceptable */
 n = 0;
 if (err[1] < err[n]) n = 1;
 if (err[2] < err[n]) n = 2;
 if (err[n] < ACCEPTABLE_EPS * Math.abs(value[n]) || err[n] < ACCEPTABLE_ATOL) {
 return value[n];
 }

 /* Maybe it really is an overflow? */
 tmp = -(0, _lgam.lgam)(v + 1.5) + (v + 1) * Math.log(z / 2);
 if (!isH) {
 tmp = Math.abs(tmp);
 }
 if (tmp > 700) {
 // sf_error("struve", SF_ERROR_OVERFLOW, "overflow in series");
 return Infinity * (0, _gammasgn.gammasgn)(v + 1.5);
 }

 /* Failure */
 // sf_error("struve", SF_ERROR_NO_RESULT, "total loss of precision");
 return NaN;
}

/*
* Power series for Struve H and L
* http://dlmf.nist.gov/11.2.1
*
* Starts to converge roughly at |n| > |z|
*/
function struvePowerSeries(v, z, isH) {
 var n = void 0,
 sgn = void 0,
 term = void 0,
 sum = void 0,
 maxterm = void 0,
 scaleexp = void 0,
 tmp = void 0,
 cterm = void 0,
 csum = void 0,
 cdiv = void 0,
 z2 = void 0,
 c2v = void 0,
 ctmp = void 0,
 err = void 0;

 if (isH) {
 sgn = -1;
 } else {
 sgn = 1;
 }

 tmp = -(0, _lgam.lgam)(v + 1.5) + (v + 1) * Math.log(z / 2);
 if (tmp < -600 || tmp > 600) {
 /* Scale exponent to postpone underflow/overflow */
 scaleexp = tmp / 2;
 tmp -= scaleexp;
 } else {
 scaleexp = 0;
 }

 term = 2 / Math.sqrt(Math.PI) * Math.exp(tmp) * (0, _gammasgn.gammasgn)(v + 1.5);
 sum = term;
 maxterm = 0;

 cterm = _dd2.default.create(term);
 csum = _dd2.default.create(sum);
 z2 = _dd2.default.create(sgn * z * z);
 c2v = _dd2.default.create(2 * v);

 for (n = 0; n < MAXITER; ++n) {
 /* cdiv = (3 + 2*n) * (3 + 2*n + 2*v)) */
 cdiv = _dd2.default.create(3 + 2 * n);
 ctmp = _dd2.default.create(3 + 2 * n);
 ctmp = _dd2.default.add(ctmp, c2v);
 cdiv = _dd2.default.mul(cdiv, ctmp);

 /* cterm *= z2 / cdiv */
 cterm = _dd2.default.mul(cterm, z2);
 cterm = _dd2.default.div(cterm, cdiv);

 csum = _dd2.default.add(csum, cterm);

 term = _dd2.default.toDouble(cterm);
 sum = _dd2.default.toDouble(csum);

 if (Math.abs(term) > maxterm) {
 maxterm = Math.abs(term);
 }
 if (Math.abs(term) < SUM_TINY * Math.abs(sum) || term === 0 || !isFinite(sum)) {
 break;
 }
 }

 err = Math.abs(term) + Math.abs(maxterm) * 1e-22;

 if (scaleexp !== 0) {
 sum *= Math.exp(scaleexp);
 err *= Math.exp(scaleexp);
 }

 if (sum === 0 && term === 0 && v < 0 && !isH) {
 /* Spurious underflow */
 err = Infinity;
 return [NaN, err];
 }

 return [sum, err];
}

/*
* Bessel series
* http://dlmf.nist.gov/11.4.19
*/
function struveBesselSeries(v, z, isH) {
 var n = void 0,
 term = void 0,
 cterm = void 0,
 sum = void 0,
 maxterm = void 0,
 err = void 0;

 if (isH && v < 0) {
 /* Works less reliably in this region */
 err = Infinity;
 return [NaN, err];
 }

 sum = 0;
 maxterm = 0;

 cterm = Math.sqrt(z / (2 * Math.PI));

 for (n = 0; n < MAXITER; ++n) {
 if (isH) {
 term = cterm * _cephes2.default.jv(n + v + 0.5, z) / (n + 0.5);
 cterm *= z / 2 / (n + 1);
 } else {
 term = cterm * _cephes2.default.iv(n + v + 0.5, z) / (n + 0.5);
 cterm *= -z / 2 / (n + 1);
 }
 sum += term;
 if (Math.abs(term) > maxterm) {
 maxterm = Math.abs(term);
 }
 if (Math.abs(term) < SUM_EPS * Math.abs(sum) || term === 0 || !isFinite(sum)) {
 break;
 }
 }

 err = Math.abs(term) + Math.abs(maxterm) * 1e-16;

 /* Account for potential underflow of the Bessel functions */
 err += 1e-300 * Math.abs(cterm);

 return [sum, err];
}

/*
* Large-z expansion for Struve H and L
* http://dlmf.nist.gov/11.6.1
*/
function struveAsympLargeZ(v, z, isH) {
 var n = void 0,
 sgn = void 0,
 maxiter = void 0,
 term = void 0,
 sum = void 0,
 maxterm = void 0,
 m = void 0,
 err = void 0;

 if (isH) {
 sgn = -1;
 } else {
 sgn = 1;
 }

 /* Asymptotic expansion divergenge point */
 m = z / 2;
 if (m <= 0) {
 maxiter = 0;
 } else if (m > MAXITER) {
 maxiter = MAXITER;
 } else {
 maxiter = Math.trunc(m);
 }
 if (maxiter === 0) {
 err = Infinity;
 return [NaN, err];
 }

 if (z < v) {
 /* Exclude regions where our error estimation fails */
 err = Infinity;
 return [NaN, err];
 }

 /* Evaluate sum */
 term = -sgn / Math.sqrt(Math.PI) * Math.exp(-(0, _lgam.lgam)(v + 0.5) + (v - 1) * Math.log(z / 2)) * (0, _gammasgn.gammasgn)(v + 0.5);
 sum = term;
 maxterm = 0;

 for (n = 0; n < maxiter; ++n) {
 term *= sgn * (1 + 2 * n) * (1 + 2 * n - 2 * v) / (z * z);
 sum += term;
 if (Math.abs(term) > maxterm) {
 maxterm = Math.abs(term);
 }
 if (Math.abs(term) < SUM_EPS * Math.abs(sum) || term === 0 || !isFinite(sum)) {
 break;
 }
 }

 if (isH) {
 sum += _cephes2.default.yv(v, z);
 } else {
 sum += _cephes2.default.iv(v, z);
 }

 /*
 * This error estimate is strictly speaking valid only for
 * n > v - 0.5, but numerical results indicate that it works
 * reasonably.
 */
 err = Math.abs(term) + Math.abs(maxterm) * 1e-16;

 return [sum, err];
}

exports.struveH = struveH;
exports.struveL = struveL;
},{"../cephes.js":48,"../cephes/dd.js":52,"../cephes/gamma.js":55,"../cephes/gamma/lgam.js":56,"./gammasgn.js":46}],48:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});

var _createClass = function () { function defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } } return function (Constructor, protoProps, staticProps) { if (protoProps) defineProperties(Constructor.prototype, protoProps); if (staticProps) defineProperties(Constructor, staticProps); return Constructor; }; }();

var _hyp2f = require('./cephes/hyp2f1.js');

var hyp2f1Cephes = _interopRequireWildcard(_hyp2f);

var _gamma = require('./cephes/gamma.js');

var gammaCephes = _interopRequireWildcard(_gamma);

var _j = require('./cephes/j0.js');

var j0Cephes = _interopRequireWildcard(_j);

var _j2 = require('./cephes/j1.js');

var j1Cephes = _interopRequireWildcard(_j2);

var _jv = require('./cephes/jv.js');

var jvCephes = _interopRequireWildcard(_jv);

var _yn = require('./cephes/yn.js');

var ynCephes = _interopRequireWildcard(_yn);

var _yv = require('./cephes/yv.js');

var yvCephes = _interopRequireWildcard(_yv);

var _i = require('./cephes/i0.js');

var i0Cephes = _interopRequireWildcard(_i);

var _i2 = require('./cephes/i1.js');

var i1Cephes = _interopRequireWildcard(_i2);

var _iv = require('./cephes/iv.js');

var ivCephes = _interopRequireWildcard(_iv);

var _k = require('./cephes/k0.js');

var k0Cephes = _interopRequireWildcard(_k);

var _k2 = require('./cephes/k1.js');

var k1Cephes = _interopRequireWildcard(_k2);

var _kn = require('./cephes/kn.js');

var knCephes = _interopRequireWildcard(_kn);

function _interopRequireWildcard(obj) { if (obj && obj.__esModule) { return obj; } else { var newObj = {}; if (obj != null) { for (var key in obj) { if (Object.prototype.hasOwnProperty.call(obj, key)) newObj[key] = obj[key]; } } newObj.default = obj; return newObj; } }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

// TODO discuss this pattern, we just want a little namespacing in the long run
// this is just a quick and dirty way to get some before we design the API
var Cephes = function () {
 function Cephes() {
 _classCallCheck(this, Cephes);
 }

 _createClass(Cephes, null, [{
 key: 'hyp2f1',
 value: function hyp2f1(a, b, c, x) {
 return hyp2f1Cephes.hyp2f1(a, b, c, x);
 }
 }, {
 key: 'gamma',
 value: function gamma(x) {
 return gammaCephes.gamma(x);
 }
 }, {
 key: 'j0',
 value: function j0(x) {
 return j0Cephes.j0(x);
 }
 }, {
 key: 'j1',
 value: function j1(x) {
 return j1Cephes.j1(x);
 }
 }, {
 key: 'jv',
 value: function jv(n, x) {
 return jvCephes.jv(n, x);
 }
 }, {
 key: 'y0',
 value: function y0(x) {
 return j0Cephes.y0(x);
 }
 }, {
 key: 'y1',
 value: function y1(x) {
 return j1Cephes.y1(x);
 }
 }, {
 key: 'yn',
 value: function yn(n, x) {
 return ynCephes.yn(n, x);
 }
 }, {
 key: 'yv',
 value: function yv(n, x) {
 return yvCephes.yv(n, x);
 }
 }, {
 key: 'i0',
 value: function i0(x) {
 return i0Cephes.i0(x);
 }
 }, {
 key: 'ie',
 value: function ie(x) {
 return i0Cephes.i0e(x);
 }
 }, {
 key: 'i1',
 value: function i1(x) {
 return i1Cephes.i1(x);
 }
 }, {
 key: 'i1e',
 value: function i1e(x) {
 return i1Cephes.i1e(x);
 }
 }, {
 key: 'iv',
 value: function iv(v, x) {
 return ivCephes.iv(v, x);
 }
 }, {
 key: 'k0',
 value: function k0(x) {
 return k0Cephes.k0(x);
 }
 }, {
 key: 'k0e',
 value: function k0e(x) {
 return k0Cephes.k0e(x);
 }
 }, {
 key: 'k1',
 value: function k1(x) {
 return k1Cephes.k1(x);
 }
 }, {
 key: 'k1e',
 value: function k1e(x) {
 return k1Cephes.k1e(x);
 }
 }, {
 key: 'kn',
 value: function kn(n, x) {
 return knCephes.kn(n, x);
 }
 }]);

 return Cephes;
}();

exports.default = Cephes;
},{"./cephes/gamma.js":55,"./cephes/hyp2f1.js":57,"./cephes/i0.js":62,"./cephes/i1.js":63,"./cephes/iv.js":64,"./cephes/j0.js":67,"./cephes/j1.js":68,"./cephes/jv.js":69,"./cephes/k0.js":75,"./cephes/k1.js":76,"./cephes/kn.js":77,"./cephes/yn.js":80,"./cephes/yv.js":81}],49:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.airy = airy;

var _polevl = require('./polevl.js');

var _constants = require('./constants.js');

var constants = _interopRequireWildcard(_constants);

function _interopRequireWildcard(obj) { if (obj && obj.__esModule) { return obj; } else { var newObj = {}; if (obj != null) { for (var key in obj) { if (Object.prototype.hasOwnProperty.call(obj, key)) newObj[key] = obj[key]; } } newObj.default = obj; return newObj; } }

// translator's note:
// return ai-bip as array + return flag to match original (-1 = error, 0 = normal)
/**
 * @file airy.js Airy function
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, ai, aip, bi, bip;
 * int airy();
 *
 * airy( x, _&ai, _&aip, _&bi, _&bip );
 *
 *
 *
 * DESCRIPTION:
 *
 * Solution of the differential equation
 *
 * y"(x) = xy.
 *
 * The function returns the two independent solutions Ai, Bi
 * and their first derivatives Ai'(x), Bi'(x).
 *
 * Evaluation is by power series summation for small x,
 * by rational minimax approximations for large x.
 *
 *
 *
 * ACCURACY:
 * Error criterion is absolute when function <= 1, relative
 * when function > 1, except * denotes relative error criterion.
 * For large negative x, the absolute error increases as x^1.5.
 * For large positive x, the relative error increases as x^1.5.
 *
 * Arithmetic domain function # trials peak rms
 * IEEE -10, 0 Ai 10000 1.6e-15 2.7e-16
 * IEEE 0, 10 Ai 10000 2.3e-14* 1.8e-15*
 * IEEE -10, 0 Ai' 10000 4.6e-15 7.6e-16
 * IEEE 0, 10 Ai' 10000 1.8e-14* 1.5e-15*
 * IEEE -10, 10 Bi 30000 4.2e-15 5.3e-16
 * IEEE -10, 10 Bi' 30000 4.9e-15 7.3e-16
 *
 *
 * Cephes Math Library Release 2.8: June, 2000
 * Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
 * Ported to ECMAScript 2018
 * Copyright (c) 2018, Kings Distributed Systems
 *
 * @author KC Erb, kc@kcerb.com
 * @date April 2018
 */
function airy(x, ai, aip, bi, bip) {
 var c1 = 0.35502805388781723926;
 var c2 = 0.258819403792806798405;
 var sqrt3 = 1.732050807568877293527;
 var sqpii = 5.64189583547756286948E-1;
 var MAXAIRY = 103.892;

 var AN = [3.46538101525629032477E-1, 1.20075952739645805542E1, 7.62796053615234516538E1, 1.68089224934630576269E2, 1.59756391350164413639E2, 7.05360906840444183113E1, 1.40264691163389668864E1, 9.99999999999999995305E-1];

 var AD = [5.67594532638770212846E-1, 1.47562562584847203173E1, 8.45138970141474626562E1, 1.77318088145400459522E2, 1.64234692871529701831E2, 7.14778400825575695274E1, 1.40959135607834029598E1, 1.00000000000000000470E0];

 var APN = [6.13759184814035759225E-1, 1.47454670787755323881E1, 8.20584123476060982430E1, 1.71184781360976385540E2, 1.59317847137141783523E2, 6.99778599330103016170E1, 1.39470856980481566958E1, 1.00000000000000000550E0];

 var APD = [3.34203677749736953049E-1, 1.11810297306158156705E1, 7.11727352147859965283E1, 1.58778084372838313640E2, 1.53206427475809220834E2, 6.86752304592780337944E1, 1.38498634758259442477E1, 9.99999999999999994502E-1];

 var BN16 = [-2.53240795869364152689E-1, 5.75285167332467384228E-1, -3.29907036873225371650E-1, 6.44404068948199951727E-2, -3.82519546641336734394E-3];

 var BD16 = [
 /* 1.00000000000000000000E0, */
 -7.15685095054035237902E0, 1.06039580715664694291E1, -5.23246636471251500874E0, 9.57395864378383833152E-1, -5.50828147163549611107E-2];

 var BPPN = [4.65461162774651610328E-1, -1.08992173800493920734E0, 6.38800117371827987759E-1, -1.26844349553102907034E-1, 7.62487844342109852105E-3];

 var BPPD = [
 /* 1.00000000000000000000E0, */
 -8.70622787633159124240E0, 1.38993162704553213172E1, -7.14116144616431159572E0, 1.34008595960680518666E0, -7.84273211323341930448E-2];

 var AFN = [-1.31696323418331795333E-1, -6.26456544431912369773E-1, -6.93158036036933542233E-1, -2.79779981545119124951E-1, -4.91900132609500318020E-2, -4.06265923594885404393E-3, -1.59276496239262096340E-4, -2.77649108155232920844E-6, -1.67787698489114633780E-8];

 var AFD = [
 /* 1.00000000000000000000E0, */
 1.33560420706553243746E1, 3.26825032795224613948E1, 2.67367040941499554804E1, 9.18707402907259625840E0, 1.47529146771666414581E0, 1.15687173795188044134E-1, 4.40291641615211203805E-3, 7.54720348287414296618E-5, 4.51850092970580378464E-7];

 var AGN = [1.97339932091685679179E-2, 3.91103029615688277255E-1, 1.06579897599595591108E0, 9.39169229816650230044E-1, 3.51465656105547619242E-1, 6.33888919628925490927E-2, 5.85804113048388458567E-3, 2.82851600836737019778E-4, 6.98793669997260967291E-6, 8.11789239554389293311E-8, 3.41551784765923618484E-10];

 var AGD = [
 /* 1.00000000000000000000E0, */
 9.30892908077441974853E0, 1.98352928718312140417E1, 1.55646628932864612953E1, 5.47686069422975497931E0, 9.54293611618961883998E-1, 8.64580826352392193095E-2, 4.12656523824222607191E-3, 1.01259085116509135510E-4, 1.17166733214413521882E-6, 4.91834570062930015649E-9];

 var APFN = [1.85365624022535566142E-1, 8.86712188052584095637E-1, 9.87391981747398547272E-1, 4.01241082318003734092E-1, 7.10304926289631174579E-2, 5.90618657995661810071E-3, 2.33051409401776799569E-4, 4.08718778289035454598E-6, 2.48379932900442457853E-8];

 var APFD = [
 /* 1.00000000000000000000E0, */
 1.47345854687502542552E1, 3.75423933435489594466E1, 3.14657751203046424330E1, 1.09969125207298778536E1, 1.78885054766999417817E0, 1.41733275753662636873E-1, 5.44066067017226003627E-3, 9.39421290654511171663E-5, 5.65978713036027009243E-7];

 var APGN = [-3.55615429033082288335E-2, -6.37311518129435504426E-1, -1.70856738884312371053E0, -1.50221872117316635393E0, -5.63606665822102676611E-1, -1.02101031120216891789E-1, -9.48396695961445269093E-3, -4.60325307486780994357E-4, -1.14300836484517375919E-5, -1.33415518685547420648E-7, -5.63803833958893494476E-10];

 var APGD = [
 /* 1.00000000000000000000E0, */
 9.85865801696130355144E0, 2.16401867356585941885E1, 1.73130776389749389525E1, 6.17872175280828766327E0, 1.08848694396321495475E0, 9.95005543440888479402E-2, 4.78468199683886610842E-3, 1.18159633322838625562E-4, 1.37480673554219441465E-6, 5.79912514929147598821E-9];

 var z = void 0,
 zz = void 0,
 t = void 0,
 f = void 0,
 g = void 0,
 uf = void 0,
 ug = void 0,
 k = void 0,
 zeta = void 0,
 theta = void 0,
 domflg = void 0;

 domflg = 0;
 if (x > MAXAIRY) {
 return [ai, aip, bi, bip, -1];
 }

 if (x < -2.09) {
 domflg = 15;
 t = Math.sqrt(-x);
 zeta = -2.0 * x * t / 3.0;
 t = Math.sqrt(t);
 k = sqpii / t;
 z = 1.0 / zeta;
 zz = z * z;
 uf = 1.0 + zz * (0, _polevl.polevl)(zz, AFN, 8) / (0, _polevl.p1evl)(zz, AFD, 9);
 ug = z * (0, _polevl.polevl)(zz, AGN, 10) / (0, _polevl.p1evl)(zz, AGD, 10);
 theta = zeta + 0.25 * Math.PI;
 f = Math.sin(theta);
 g = Math.cos(theta);
 ai = k * (f * uf - g * ug);
 bi = k * (g * uf + f * ug);
 uf = 1.0 + zz * (0, _polevl.polevl)(zz, APFN, 8) / (0, _polevl.p1evl)(zz, APFD, 9);
 ug = z * (0, _polevl.polevl)(zz, APGN, 10) / (0, _polevl.p1evl)(zz, APGD, 10);
 k = sqpii * t;
 aip = -k * (g * uf + f * ug);
 bip = k * (f * uf - g * ug);
 return [ai, aip, bi, bip, 0];
 }

 if (x >= 2.09) {
 /* cbrt(9) */
 domflg = 5;
 t = Math.sqrt(x);
 zeta = 2.0 * x * t / 3.0;
 g = Math.exp(zeta);
 t = Math.sqrt(t);
 k = 2.0 * t * g;
 z = 1.0 / zeta;
 f = (0, _polevl.polevl)(z, AN, 7) / (0, _polevl.polevl)(z, AD, 7);
 ai = sqpii * f / k;
 k = -0.5 * sqpii * t / g;
 f = (0, _polevl.polevl)(z, APN, 7) / (0, _polevl.polevl)(z, APD, 7);
 aip = f * k;

 if (x > 8.3203353) {
 /* zeta > 16 */
 f = z * (0, _polevl.polevl)(z, BN16, 4) / (0, _polevl.p1evl)(z, BD16, 5);
 k = sqpii * g;
 bi = k * (1.0 + f) / t;
 f = z * (0, _polevl.polevl)(z, BPPN, 4) / (0, _polevl.p1evl)(z, BPPD, 5);
 bip = k * t * (1.0 + f);
 return [ai, aip, bi, bip, 0];
 }
 }

 f = 1.0;
 g = x;
 t = 1.0;
 uf = 1.0;
 ug = x;
 k = 1.0;
 z = x * x * x;
 while (t > constants.MACHEP) {
 uf *= z;
 k += 1.0;
 uf /= k;
 ug *= z;
 k += 1.0;
 ug /= k;
 uf /= k;
 f += uf;
 k += 1.0;
 ug /= k;
 g += ug;
 t = Math.abs(uf / f);
 }
 uf = c1 * f;
 ug = c2 * g;

 if ((domflg & 1) === 0) ai = uf - ug;
 if ((domflg & 2) === 0) bi = sqrt3 * (uf + ug);

 /* the deriviative of ai */
 k = 4.0;
 uf = x * x / 2.0;
 ug = z / 3.0;
 f = uf;
 g = 1.0 + ug;
 uf /= 3.0;
 t = 1.0;

 while (t > constants.MACHEP) {
 uf *= z;
 ug /= k;
 k += 1.0;
 ug *= z;
 uf /= k;
 f += uf;
 k += 1.0;
 ug /= k;
 uf /= k;
 g += ug;
 k += 1.0;
 t = Math.abs(ug / g);
 }

 uf = c1 * f;
 ug = c2 * g;
 if ((domflg & 4) === 0) aip = uf - ug;
 if ((domflg & 8) === 0) bip = sqrt3 * (uf + ug);
 return [ai, aip, bi, bip, 0];
}
},{"./constants.js":51,"./polevl.js":78}],50:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
/**
 * @file chbevl.js Evaluate Chebyshev series
 *
 *
 *
 *
 * SYNOPSIS:
 *
 * int N;
 * double x, y, coef[N], chebevl();
 *
 * y = chbevl( x, coef, N );
 *
 *
 *
 * DESCRIPTION:
 *
 * Evaluates the series
 *
 * N-1
 * - '
 * y = > coef[i] T (x/2)
 * - i
 * i=0
 *
 * of Chebyshev polynomials Ti at argument x/2.
 *
 * Coefficients are stored in reverse order, i.e. the zero
 * order term is last in the array. Note N is the number of
 * coefficients, not the order.
 *
 * If coefficients are for the interval a to b, x must
 * have been transformed to x -> 2(2x - b - a)/(b-a) before
 * entering the routine. This maps x from (a, b) to (-1, 1),
 * over which the Chebyshev polynomials are defined.
 *
 * If the coefficients are for the inverted interval, in
 * which (a, b) is mapped to (1/b, 1/a), the transformation
 * required is x -> 2(2ab/x - b - a)/(b-a). If b is infinity,
 * this becomes x -> 4a/x - 1.
 *
 *
 *
 * SPEED:
 *
 * Taking advantage of the recurrence properties of the
 * Chebyshev polynomials, the routine requires one more
 * addition per loop than evaluating a nested polynomial of
 * the same degree.
 *
 *
 * Cephes Math Library Release 2.0: April, 1987
 * Copyright 1985, 1987 by Stephen L. Moshier
 * Ported to ECMAScript 2018
 * Copyright (c) 2018, Kings Distributed Systems
 *
 * @author KC Erb, kc@kcerb.com
 * @date April 2018
 */
function chbevl(x, array, n) {
 var b0 = void 0,
 b1 = void 0,
 b2 = void 0,
 i = void 0;

 b0 = array[0];
 b1 = 0.0;

 for (i = 1; i < n; i++) {
 b2 = b1;
 b1 = b0;
 b0 = x * b1 - b2 + array[i];
 }

 return 0.5 * (b0 - b2);
}

exports.chbevl = chbevl;
},{}],51:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
var MACHEP = exports.MACHEP = Number.EPSILON / 2;
var DBL_MAX = exports.DBL_MAX = Number.MAX_VALUE;
var MAXLOG = exports.MAXLOG = 7.09782712893383996843E2;
var MINLOG = exports.MINLOG = -7.08396418532264106224E2;
var SQ2OPI = exports.SQ2OPI = 7.9788456080286535587989E-1; // sqrt( 2/pi )
var LOGSQ2 = exports.LOGSQ2 = 3.46573590279972654709E-1; // log(2)/2
var THPIO4 = exports.THPIO4 = 2.35619449019234492885; // 3*pi/4
var NPY_PI_4 = exports.NPY_PI_4 = 0.78539816339744830962; // pi / 4
var NPY_2_PI = exports.NPY_2_PI = 0.63661977236758134308; // 2 / pi
var MAXGAM = exports.MAXGAM = 171.624376956302725;
var SQTPI = exports.SQTPI = 2.50662827463100050242E0;
var LS2PI = exports.LS2PI = 0.918938533204672; /* log( sqrt( 2*pi ) ) */
var EULER = exports.EULER = 0.57721566490153286061;
var MAXITER = exports.MAXITER = 500;
},{}],52:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});

var _createClass = function () { function defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } } return function (Constructor, protoProps, staticProps) { if (protoProps) defineProperties(Constructor.prototype, protoProps); if (staticProps) defineProperties(Constructor, staticProps); return Constructor; }; }(); // This is a (currently) brief translation of the Double-double precision
// arithmetic package in cephes. For now we'll add the methods we need and let
// it grow naturally.


var _ddRealIdefs = require('./dd/ddRealIdefs.js');

var ddRealIdefs = _interopRequireWildcard(_ddRealIdefs);

function _interopRequireWildcard(obj) { if (obj && obj.__esModule) { return obj; } else { var newObj = {}; if (obj != null) { for (var key in obj) { if (Object.prototype.hasOwnProperty.call(obj, key)) newObj[key] = obj[key]; } } newObj.default = obj; return newObj; } }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

var DD = function () {
 function DD() {
 _classCallCheck(this, DD);
 }

 _createClass(DD, null, [{
 key: 'create',
 value: function create(hi) {
 var low = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : 0;
 return [hi, low];
 }
 }, {
 key: 'toDouble',
 value: function toDouble(a) {
 return a[0];
 }
 // TODO - let these accept non double-double args or at least do error checking
 // that they are both double doubles (i.e. arrays)

 }, {
 key: 'add',
 value: function add(a, b) {
 return ddRealIdefs.ddAdd(a, b);
 }
 }, {
 key: 'sub',
 value: function sub(a, b) {
 return ddRealIdefs.ddSub(a, b);
 }
 }, {
 key: 'mul',
 value: function mul(a, b) {
 return ddRealIdefs.ddMul(a, b);
 }
 }, {
 key: 'div',
 value: function div(a, b) {
 return ddRealIdefs.ddDiv(a, b);
 }
 }]);

 return DD;
}();

exports.default = DD;
},{"./dd/ddRealIdefs.js":54}],53:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});

var _slicedToArray = function () { function sliceIterator(arr, i) { var _arr = []; var _n = true; var _d = false; var _e = undefined; try { for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) { _arr.push(_s.value); if (i && _arr.length === i) break; } } catch (err) { _d = true; _e = err; } finally { try { if (!_n && _i["return"]) _i["return"](); } finally { if (_d) throw _e; } } return _arr; } return function (arr, i) { if (Array.isArray(arr)) { return arr; } else if (Symbol.iterator in Object(arr)) { return sliceIterator(arr, i); } else { throw new TypeError("Invalid attempt to destructure non-iterable instance"); } }; }();

/*
 * Ported to ECMAScript 2018 by KC Erb of Kings Distributed Systems from
 * original C code which contained the following:
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract numbers DE-AC03-76SF00098 and
 * DE-AC02-05CH11231.
 *
 * Copyright (c) 2003-2009, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from U.S. Dept. of Energy) All rights reserved.
 *
 * By downloading or using this software you are agreeing to the modified
 * BSD license "BSD-LBNL-License.doc" (see LICENSE.txt).
 */
/*
 * Contains small functions (suitable for inlining) in the double-double
 * arithmetic package.
 */
var _DD_SPLIT_THRESH = 6.69692879491417e+299; // = 2^996
var _DD_SPLITTER = 134217729.0; // = 2^27 + 1

function twoSum(a, b) {
 var s = a + b;
 var c = s - a;
 var d = b - c;
 var e = s - c;
 var err = a - e + d;
 return [s, err];
}

function quickTwoSum(a, b) {
 var s = a + b;
 var c = s - a;
 var err = b - c;
 return [s, err];
}

/* Computes fl(a*b) and err(a*b). */
function twoProd(a, b) {
 var aHi = void 0,
 aLo = void 0,
 bHi = void 0,
 bLo = void 0,
 c = void 0,
 d = void 0,
 err = void 0;
 var p = a * b;

 var _twoSplit = twoSplit(a);

 var _twoSplit2 = _slicedToArray(_twoSplit, 2);

 aHi = _twoSplit2[0];
 aLo = _twoSplit2[1];

 var _twoSplit3 = twoSplit(b);

 var _twoSplit4 = _slicedToArray(_twoSplit3, 2);

 bHi = _twoSplit4[0];
 bLo = _twoSplit4[1];

 c = aHi * bHi - p;
 d = c + aHi * bLo + aLo * bHi;
 err = d + aLo * bLo;
 return [p, err];
}

/* Computes high word and lo word of a */
function twoSplit(a) {
 var temp = void 0,
 tempma = void 0,
 hi = void 0,
 lo = void 0;
 if (a > _DD_SPLIT_THRESH || a < -_DD_SPLIT_THRESH) {
 a *= 3.7252902984619140625e-09; // 2^-28
 temp = _DD_SPLITTER * a;
 tempma = temp - a;
 hi = temp - tempma;
 lo = a - hi;
 hi *= 268435456.0; // 2^28
 lo *= 268435456.0; // 2^28
 } else {
 temp = _DD_SPLITTER * a;
 tempma = temp - a;
 hi = temp - tempma;
 lo = a - hi;
 }
 return [hi, lo];
}

exports.twoSum = twoSum;
exports.quickTwoSum = quickTwoSum;
exports.twoProd = twoProd;
},{}],54:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.ddDiv = exports.ddMul = exports.ddSub = exports.ddAdd = undefined;

var _slicedToArray = function () { function sliceIterator(arr, i) { var _arr = []; var _n = true; var _d = false; var _e = undefined; try { for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) { _arr.push(_s.value); if (i && _arr.length === i) break; } } catch (err) { _d = true; _e = err; } finally { try { if (!_n && _i["return"]) _i["return"](); } finally { if (_d) throw _e; } } return _arr; } return function (arr, i) { if (Array.isArray(arr)) { return arr; } else if (Symbol.iterator in Object(arr)) { return sliceIterator(arr, i); } else { throw new TypeError("Invalid attempt to destructure non-iterable instance"); } }; }(); /*
 * Ported to ECMAScript 2018 by KC Erb of Kings Distributed Systems from
 * original C code which contained the following:
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract numbers DE-AC03-76SF00098 and
 * DE-AC02-05CH11231.
 *
 * Copyright (c) 2003-2009, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from U.S. Dept. of Energy) All rights reserved.
 *
 * By downloading or using this software you are agreeing to the modified
 * BSD license "BSD-LBNL-License.doc" (see LICENSE.txt).
 */
/*
 * Contains small functions (suitable for inlining) in the double-double
 * arithmetic package.
 */


var _ddIdefs = require('./ddIdefs.js');

var ddIdefs = _interopRequireWildcard(_ddIdefs);

function _interopRequireWildcard(obj) { if (obj && obj.__esModule) { return obj; } else { var newObj = {}; if (obj != null) { for (var key in obj) { if (Object.prototype.hasOwnProperty.call(obj, key)) newObj[key] = obj[key]; } } newObj.default = obj; return newObj; } }

function ddNeg(a) {
 return [-a[0], -a[1]];
}

// dd_ieee_add in orig
function ddAdd(a, b) {
 var s1 = void 0,
 s2 = void 0,
 t1 = void 0,
 t2 = void 0;

 var _ddIdefs$twoSum = ddIdefs.twoSum(a[0], b[0]);

 var _ddIdefs$twoSum2 = _slicedToArray(_ddIdefs$twoSum, 2);

 s1 = _ddIdefs$twoSum2[0];
 s2 = _ddIdefs$twoSum2[1];

 var _ddIdefs$twoSum3 = ddIdefs.twoSum(a[1], b[1]);

 var _ddIdefs$twoSum4 = _slicedToArray(_ddIdefs$twoSum3, 2);

 t1 = _ddIdefs$twoSum4[0];
 t2 = _ddIdefs$twoSum4[1];

 s2 += t1;

 var _ddIdefs$quickTwoSum = ddIdefs.quickTwoSum(s1, s2);

 var _ddIdefs$quickTwoSum2 = _slicedToArray(_ddIdefs$quickTwoSum, 2);

 s1 = _ddIdefs$quickTwoSum2[0];
 s2 = _ddIdefs$quickTwoSum2[1];

 s2 += t2;

 var _ddIdefs$quickTwoSum3 = ddIdefs.quickTwoSum(s1, s2);

 var _ddIdefs$quickTwoSum4 = _slicedToArray(_ddIdefs$quickTwoSum3, 2);

 s1 = _ddIdefs$quickTwoSum4[0];
 s2 = _ddIdefs$quickTwoSum4[1];

 return [s1, s2];
}

// first arg double2[], second arg double
function ddAddDdD(a, b) {
 var s1 = void 0,
 s2 = void 0;

 var _ddIdefs$twoSum5 = ddIdefs.twoSum(a[0], b);

 var _ddIdefs$twoSum6 = _slicedToArray(_ddIdefs$twoSum5, 2);

 s1 = _ddIdefs$twoSum6[0];
 s2 = _ddIdefs$twoSum6[1];

 s2 += a[1];

 var _ddIdefs$quickTwoSum5 = ddIdefs.quickTwoSum(s1, s2);

 var _ddIdefs$quickTwoSum6 = _slicedToArray(_ddIdefs$quickTwoSum5, 2);

 s1 = _ddIdefs$quickTwoSum6[0];
 s2 = _ddIdefs$quickTwoSum6[1];

 return [s1, s2];
}

function ddSub(a, b) {
 return ddAdd(a, ddNeg(b));
}

function ddMul(a, b) {
 var p1 = void 0,
 p2 = void 0;

 var _ddIdefs$twoProd = ddIdefs.twoProd(a[0], b[0]);

 var _ddIdefs$twoProd2 = _slicedToArray(_ddIdefs$twoProd, 2);

 p1 = _ddIdefs$twoProd2[0];
 p2 = _ddIdefs$twoProd2[1];

 p2 += a[0] * b[1] + a[1] * b[0];

 var _ddIdefs$quickTwoSum7 = ddIdefs.quickTwoSum(p1, p2);

 var _ddIdefs$quickTwoSum8 = _slicedToArray(_ddIdefs$quickTwoSum7, 2);

 p1 = _ddIdefs$quickTwoSum8[0];
 p2 = _ddIdefs$quickTwoSum8[1];

 return [p1, p2];
}

// first arg double2[], second arg double
function ddMulDdD(a, b) {
 var p1 = void 0,
 p2 = void 0,
 e1 = void 0,
 e2 = void 0;

 var _ddIdefs$twoProd3 = ddIdefs.twoProd(a[0], b);

 var _ddIdefs$twoProd4 = _slicedToArray(_ddIdefs$twoProd3, 2);

 p1 = _ddIdefs$twoProd4[0];
 e1 = _ddIdefs$twoProd4[1];

 var _ddIdefs$twoProd5 = ddIdefs.twoProd(a[1], b);

 var _ddIdefs$twoProd6 = _slicedToArray(_ddIdefs$twoProd5, 2);

 p2 = _ddIdefs$twoProd6[0];
 e2 = _ddIdefs$twoProd6[1];

 var _ddIdefs$quickTwoSum9 = ddIdefs.quickTwoSum(p1, e2 + p2 + e1);

 var _ddIdefs$quickTwoSum10 = _slicedToArray(_ddIdefs$quickTwoSum9, 2);

 p1 = _ddIdefs$quickTwoSum10[0];
 e1 = _ddIdefs$quickTwoSum10[1];

 return [p1, e1];
}
// dd_accurate_div in original
function ddDiv(a, b) {
 var q1 = void 0,
 q2 = void 0,
 q3 = void 0,
 r = void 0;
 q1 = a[0] / b[0]; /* approximate quotient */

 r = ddSub(a, ddMulDdD(b, q1));

 q2 = r[0] / b[0];
 r = ddSub(r, ddMulDdD(b, q2));

 q3 = r[0] / b[0];

 var _ddIdefs$quickTwoSum11 = ddIdefs.quickTwoSum(q1, q2);

 var _ddIdefs$quickTwoSum12 = _slicedToArray(_ddIdefs$quickTwoSum11, 2);

 q1 = _ddIdefs$quickTwoSum12[0];
 q2 = _ddIdefs$quickTwoSum12[1];

 r = ddAddDdD([q1, q2], q3);
 return r;
}

exports.ddAdd = ddAdd;
exports.ddSub = ddSub;
exports.ddMul = ddMul;
exports.ddDiv = ddDiv;
},{"./ddIdefs.js":53}],55:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.gamma = undefined;

var _polevl = require('./polevl.js');

var _constants = require('./constants.js');

var constants = _interopRequireWildcard(_constants);

function _interopRequireWildcard(obj) { if (obj && obj.__esModule) { return obj; } else { var newObj = {}; if (obj != null) { for (var key in obj) { if (Object.prototype.hasOwnProperty.call(obj, key)) newObj[key] = obj[key]; } } newObj.default = obj; return newObj; } }

/**
 * @file gamma.js Gamma function
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, Gamma();
 *
 * y = Gamma( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Gamma function of the argument. The result is
 * correctly signed.
 *
 * Arguments |x| <= 34 are reduced by recurrence and the function
 * approximated by a rational function of degree 6/7 in the
 * interval (2,3). Large arguments are handled by Stirling's
 * formula. Large negative arguments are made positive using
 * a reflection formula.
 *
 *
 * ACCURACY:
 *
 * Relative error:
 * arithmetic domain # trials peak rms
 * IEEE -170,-33 20000 2.3e-15 3.3e-16
 * IEEE -33, 33 20000 9.4e-16 2.2e-16
 * IEEE 33, 171.6 20000 2.3e-15 3.2e-16
 *
 * Error for arguments outside the test range will be larger
 * owing to error amplification by the exponential function.
 *
 * Cephes Math Library Release 2.2: July, 1992
 * Copyright 1984, 1987, 1989, 1992 by Stephen L. Moshier
 * Ported to ECMAScript 2018
 * Copyright (c) 2018, Kings Distributed Systems
 *
 * @author KC Erb, kc@kcerb.com
 * @date April 2018
 */
var P = new Float64Array([1.60119522476751861407E-4, 1.19135147006586384913E-3, 1.04213797561761569935E-2, 4.76367800457137231464E-2, 2.07448227648435975150E-1, 4.94214826801497100753E-1, 9.99999999999999996796E-1]);

var Q = new Float64Array([-2.31581873324120129819E-5, 5.39605580493303397842E-4, -4.45641913851797240494E-3, 1.18139785222060435552E-2, 3.58236398605498653373E-2, -2.34591795718243348568E-1, 7.14304917030273074085E-2, 1.00000000000000000320E0]);

function gamma(x) {
 var p = void 0,
 q = void 0,
 z = void 0,
 i = void 0;
 var sgngam = 1;

 if (!isFinite(x)) {
 return x;
 }
 q = Math.abs(x);
 if (q > 33.0) {
 if (x < 0.0) {
 p = Math.floor(q);
 if (p === q) {
 // mtherr("Gamma", OVERFLOW);
 return Infinity;
 }
 i = p;
 if ((i & 1) === 0) {
 sgngam = -1;
 }

 z = q - p;
 if (z > 0.5) {
 p += 1.0;
 z = q - p;
 }
 z = q * Math.sin(Math.PI * z);
 if (z === 0.0) {
 return sgngam * Infinity;
 }
 z = Math.abs(z);
 z = Math.PI / (z * stirf(q));
 } else {
 z = stirf(x);
 }
 return sgngam * z;
 }

 z = 1.0;
 while (x >= 3.0) {
 x -= 1.0;
 z *= x;
 }

 while (x < 0.0) {
 if (x > -1.E-9) return small(x, z);
 z /= x;
 x += 1.0;
 }

 while (x < 2.0) {
 if (x < 1.e-9) {
 return small(x, z);
 }
 z /= x;
 x += 1.0;
 }

 if (x === 2.0) return z;

 x -= 2.0;
 p = (0, _polevl.polevl)(x, P, 6);
 q = (0, _polevl.polevl)(x, Q, 7);
 return z * p / q;
};

function small(x, z) {
 if (x === 0.0) {
 // mtherr("Gamma", OVERFLOW);
 return Infinity;
 } else {
 return z / ((1.0 + 0.5772156649015329 * x) * x);
 }
}

// Stirling's formula for the Gamma function
var STIR = [8.33333333333482257126E-2, 3.47222221605458667310E-3, -2.68132617805781232825E-3, -2.29549961613378126380E-4, 7.87311395793093628397E-4];

function stirf(x) {
 var y = void 0,
 w = void 0,
 v = void 0;

 var MAXSTIR = 143.01608;

 if (x >= constants.MAXGAM) {
 return Infinity;
 }
 w = 1.0 / x;
 w = 1.0 + w * (0, _polevl.polevl)(w, STIR, 4);
 y = Math.exp(x);

 // Avoid overflow in pow()
 if (x > MAXSTIR) {
 v = Math.pow(x, 0.5 * x - 0.25);
 y = v * (v / y);
 } else {
 y = Math.pow(x, x - 0.5) / y;
 }
 y = constants.SQTPI * y * w;
 return y;
}

exports.gamma = gamma;
},{"./constants.js":51,"./polevl.js":78}],56:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.lgamSgn = exports.lgam = undefined;

var _slicedToArray = function () { function sliceIterator(arr, i) { var _arr = []; var _n = true; var _d = false; var _e = undefined; try { for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) { _arr.push(_s.value); if (i && _arr.length === i) break; } } catch (err) { _d = true; _e = err; } finally { try { if (!_n && _i["return"]) _i["return"](); } finally { if (_d) throw _e; } } return _arr; } return function (arr, i) { if (Array.isArray(arr)) { return arr; } else if (Symbol.iterator in Object(arr)) { return sliceIterator(arr, i); } else { throw new TypeError("Invalid attempt to destructure non-iterable instance"); } }; }(); /**
 * @file lgam.js Natural logarithm of Gamma function
 *
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, lgam();
 *
 * y = lgam( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the base e (2.718...) logarithm of the absolute
 * value of the Gamma function of the argument.
 *
 * For arguments greater than 13, the logarithm of the Gamma
 * function is approximated by the logarithmic version of
 * Stirling's formula using a polynomial approximation of
 * degree 4. Arguments between -33 and +33 are reduced by
 * recurrence to the interval [2,3] of a rational approximation.
 * The cosecant reflection formula is employed for arguments
 * less than -33.
 *
 * Arguments greater than MAXLGM return NPY_INFINITY and an error
 * message. MAXLGM = 2.556348e305 for IEEE arithmetic.
 *
 *
 *
 * ACCURACY:
 *
 *
 * arithmetic domain # trials peak rms
 * IEEE 0, 3 28000 5.4e-16 1.1e-16
 * IEEE 2.718, 2.556e305 40000 3.5e-16 8.3e-17
 * The error criterion was relative when the function magnitude
 * was greater than one but absolute when it was less than one.
 *
 * The following test used the relative error criterion, though
 * at certain points the relative error could be much higher than
 * indicated.
 * IEEE -200, -4 10000 4.8e-16 1.3e-16
 *
 * Cephes Math Library Release 2.2: July, 1992
 * Copyright 1984, 1987, 1989, 1992 by Stephen L. Moshier
 * Ported to ECMAScript 2018
 * Copyright (c) 2018, Kings Distributed Systems
 *
 * @author KC Erb, kc@kcerb.com
 * @date April 2018
 */


var _polevl = require('../polevl.js');

var _constants = require('../constants.js');

var constants = _interopRequireWildcard(_constants);

function _interopRequireWildcard(obj) { if (obj && obj.__esModule) { return obj; } else { var newObj = {}; if (obj != null) { for (var key in obj) { if (Object.prototype.hasOwnProperty.call(obj, key)) newObj[key] = obj[key]; } } newObj.default = obj; return newObj; } }

/* A[]: Stirling's formula expansion of log Gamma
* B[], C[]: log Gamma function between 2 and 3
*/
var A = [8.11614167470508450300E-4, -5.95061904284301438324E-4, 7.93650340457716943945E-4, -2.77777777730099687205E-3, 8.33333333333331927722E-2];

var B = [-1.37825152569120859100E3, -3.88016315134637840924E4, -3.31612992738871184744E5, -1.16237097492762307383E6, -1.72173700820839662146E6, -8.53555664245765465627E5];

var C = [
/* 1.00000000000000000000E0, */
-3.51815701436523470549E2, -1.70642106651881159223E4, -2.20528590553854454839E5, -1.13933444367982507207E6, -2.53252307177582951285E6, -2.01889141433532773231E6];

var MAXLGM = 2.556348e305;
var LOGPI = 1.14472988584940017414;

function lgam(x) {
 var _lgamSgn = lgamSgn(x),
 _lgamSgn2 = _slicedToArray(_lgamSgn, 1),
 res = _lgamSgn2[0];

 return res;
}

// Translator's note:
// returns array of [val, sign] to replace pointer functionality from original
function lgamSgn(x) {
 var p = void 0,
 q = void 0,
 u = void 0,
 w = void 0,
 z = void 0,
 i = void 0,
 sign = void 0;
 sign = 1;

 if (!isFinite(x)) return [x, sign];

 if (x < -34.0) {
 q = -x;

 var _lgamSgn3 = lgamSgn(q, sign);

 var _lgamSgn4 = _slicedToArray(_lgamSgn3, 2);

 w = _lgamSgn4[0];
 sign = _lgamSgn4[1];

 p = Math.floor(q);
 if (p === q) {
 // mtherr("lgam", SING);
 return [Infinity, sign];
 }
 i = p;
 if ((i & 1) === 0) sign = -1;else sign = 1;
 z = q - p;
 if (z > 0.5) {
 p += 1.0;
 z = p - q;
 }
 z = q * Math.sin(Math.PI * z);
 if (z === 0.0) {
 // mtherr("lgam", SING);
 return [Infinity, sign];
 }
 /* z = log(NPY_PI) - log( z ) - w; */
 z = LOGPI - Math.log(z) - w;
 return [z, sign];
 }

 if (x < 13.0) {
 z = 1.0;
 p = 0.0;
 u = x;
 while (u >= 3.0) {
 p -= 1.0;
 u = x + p;
 z *= u;
 }
 while (u < 2.0) {
 if (u === 0.0) {
 // mtherr("lgam", SING);
 return [Infinity, sign];
 }
 z /= u;
 p += 1.0;
 u = x + p;
 }
 if (z < 0.0) {
 sign = -1;
 z = -z;
 } else {
 sign = 1;
 }
 if (u === 2.0) return [Math.log(z), sign];
 p -= 2.0;
 x = x + p;
 p = x * (0, _polevl.polevl)(x, B, 5) / (0, _polevl.p1evl)(x, C, 6);
 return [Math.log(z) + p, sign];
 }

 if (x > MAXLGM) {
 return [sign * Infinity, sign];
 }

 q = (x - 0.5) * Math.log(x) - x + constants.LS2PI;
 if (x > 1.0e8) return [q, sign];

 p = 1.0 / (x * x);
 if (x >= 1000.0) {
 q += ((7.9365079365079365079365e-4 * p - 2.7777777777777777777778e-3) * p + 0.0833333333333333333333) / x;
 } else {
 q += (0, _polevl.polevl)(p, A, 4) / x;
 }
 return [q, sign];
}

exports.lgam = lgam;
exports.lgamSgn = lgamSgn;
},{"../constants.js":51,"../polevl.js":78}],57:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.hyp2f1 = undefined;

var _constants = require('./hyp2f1/constants.js');

var constants = _interopRequireWildcard(_constants);

var _gamma = require('./gamma.js');

var _hyt2f = require('./hyp2f1/hyt2f1.js');

var _hys2f = require('./hyp2f1/hys2f1.js');

function _interopRequireWildcard(obj) { if (obj && obj.__esModule) { return obj; } else { var newObj = {}; if (obj != null) { for (var key in obj) { if (Object.prototype.hasOwnProperty.call(obj, key)) newObj[key] = obj[key]; } } newObj.default = obj; return newObj; } }

//
/**
 * @file hyp2f1.js Gauss hypergeometric function F
 * 2 1
 *
 *
 * SYNOPSIS:
 *
 * double a, b, c, x, y, hyp2f1();
 *
 * y = hyp2f1( a, b, c, x );
 *
 *
 * DESCRIPTION:
 *
 *
 * hyp2f1( a, b, c, x ) = F ( a, b; c; x )
 * 2 1
 *
 * inf.
 * - a(a+1)...(a+k) b(b+1)...(b+k) k+1
 * = 1 + > ----------------------------- x .
 * - c(c+1)...(c+k) (k+1)!
 * k = 0
 *
 * Cases addressed are
 * Tests and escapes for negative integer a, b, or c
 * Linear transformation if c - a or c - b negative integer
 * Special case c = a or c = b
 * Linear transformation for x near +1
 * Transformation for x < -0.5
 * Psi function expansion if x > 0.5 and c - a - b integer
 * Conditionally, a recurrence on c to make c-a-b > 0
 *
 * x < -1 AMS 15.3.7 transformation applied (Travis Oliphant)
 * valid for b,a,c,(b-a) !== integer and (c-a),(c-b) !== negative integer
 *
 * x >= 1 is rejected (unless special cases are present)
 *
 * The parameters a, b, c are considered to be integer
 * valued if they are within 1.0e-14 of the nearest integer
 * (1.0e-13 for IEEE arithmetic).
 *
 * ACCURACY:
 *
 *
 * Relative error (-1 < x < 1):
 * arithmetic domain # trials peak rms
 * IEEE -1,7 230000 1.2e-11 5.2e-14
 *
 * Several special cases also tested with a, b, c in
 * the range -7 to 7.
 *
 * ERROR MESSAGES:
 *
 * A "partial loss of precision" message is printed if
 * the internally estimated relative error exceeds 1^-12.
 * A "singularity" message is printed on overflow or
 * in cases not addressed (such as x < -1).
 * Cephes Math Library Release 2.8: June, 2000
 * Copyright 1984, 1987, 1992, 2000 by Stephen L. Moshier
 * Ported to ECMAScript 2018
 * Copyright (c) 2018, Kings Distributed Systems
 *
 * @author KC Erb, kc@kcerb.com
 * @date April 2018
 */
function hyp2f1(a, b, c, x) {
 var cInt = void 0;
 if (x === 0.0) {
 return 1.0;
 }

 if ((a === 0 || b === 0) && c !== 0) {
 return 1.0;
 }

 var negIntA = 0;
 var negIntB = 0;
 var negIntCaOrCb = 0;

 var errObj = { val: 0.0 };
 var xAbs = Math.abs(x);
 var s = 1.0 - x;
 var aInt = Math.round(a); /* nearest integer to a */
 var bInt = Math.round(b);

 var d = c - a - b;
 var dInt = Math.round(d);
 var p = void 0,
 q = void 0,
 y = void 0;

 // a is a negative integer
 if (a <= 0 && Math.abs(a - aInt) < constants.EPS) {
 negIntA = 1;
 }

 // b is a negative integer
 if (b <= 0 && Math.abs(b - bInt) < constants.EPS) {
 negIntB = 1;
 }

 if (d <= -1 && !(Math.abs(d - dInt) > constants.EPS && s < 0) && !(negIntA || negIntB)) {
 return Math.pow(s, d) * hyp2f1(c - a, c - b, c, x);
 }
 if (d <= 0 && x === 1 && !(negIntA || negIntB)) {
 return hypdiv();
 }

 /* 2F1(a,b;b;x) = (1-x)**(-a) */
 if (xAbs < 1.0 || x === -1.0) {
 // b = c
 if (Math.abs(b - c) < constants.EPS) {
 y = s ** -a;
 return hypdon(y, errObj);
 }

 // a = c
 if (Math.abs(a - c) < constants.EPS) {
 y = s ** -b;
 return hypdon(y, errObj);
 }
 }

 if (c <= 0.0) {
 cInt = Math.round(c);
 // c is a negative integer
 if (Math.abs(c - cInt) < constants.EPS) {
 // check if termination before explosion
 if (negIntA && aInt > cInt) {
 return hypok(a, b, c, x, errObj);
 }
 if (negIntB && bInt > cInt) {
 return hypok(a, b, c, x, errObj);
 }
 return hypdiv();
 }
 }

 // function is a polynomial
 if (negIntA || negIntB) {
 return hypok(a, b, c, x, errObj);
 }

 var t1 = Math.abs(b - a);

 if (x < -2.0 && Math.abs(t1 - Math.round(t1)) > constants.EPS) {
 // This transform has a pole for b-a integer, and
 // may produce large cancellation errors for |1/x| close to 1
 var _p = hyp2f1(a, 1 - c + a, 1 - b + a, 1.0 / x);
 var _q = hyp2f1(b, 1 - c + b, 1 - a + b, 1.0 / x);
 _p *= Math.pow(-x, -a);
 _q *= Math.pow(-x, -b);
 t1 = (0, _gamma.gamma)(c);
 s = t1 * (0, _gamma.gamma)(b - a) / ((0, _gamma.gamma)(b) * (0, _gamma.gamma)(c - a));
 y = t1 * (0, _gamma.gamma)(a - b) / ((0, _gamma.gamma)(a) * (0, _gamma.gamma)(c - b));
 return s * _p + y * _q;
 } else if (x < -1.0) {
 if (Math.abs(a) < Math.abs(b)) {
 return Math.pow(s, -a) * hyp2f1(a, c - b, c, x / (x - 1));
 } else {
 return Math.pow(s, -b) * hyp2f1(b, c - a, c, x / (x - 1));
 }
 }

 // series diverges
 if (xAbs > 1.0) {
 return hypdiv();
 }

 p = c - a;
 // nearest integer to c-a
 aInt = Math.round(p);

 // negative int c - a
 if (aInt <= 0.0 && Math.abs(p - aInt) < constants.EPS) {
 negIntCaOrCb = 1;
 }

 var r = c - b;
 // nearest integer to c-b
 bInt = Math.round(r);

 if (bInt <= 0.0 && Math.abs(r - bInt) < constants.EPS) {
 negIntCaOrCb = 1;
 }

 // nearest integer to d
 dInt = Math.round(d);
 q = Math.abs(d - dInt);

 // |x|===1.0
 if (Math.abs(xAbs - 1.0) < constants.EPS) {
 if (x > 0.0) {
 if (negIntCaOrCb) {
 return d >= 0.0 ? hypf(a, b, c, d, s, x, errObj) : hypdiv();
 }
 if (d <= 0.0) {
 hypdiv();
 }
 y = (0, _gamma.gamma)(c) * (0, _gamma.gamma)(d) / ((0, _gamma.gamma)(p) * (0, _gamma.gamma)(r));
 return hypdon(y, errObj);
 }

 if (d <= -1.0) {
 return hypdiv();
 }
 }

 // Conditionally make d > 0 by recurrence on c
 if (d < 0.0) {
 // Try the power series first
 y = (0, _hyt2f.hyt2f1)(a, b, c, x, errObj);

 if (errObj.val < constants.ETHRESH) {
 return hypdon(y, errObj);
 }

 /* Apply the recurrence if power series fails */
 errObj.val = 0.0;
 var aid = 2 - dInt;
 var e = c + aid;
 var d2 = hyp2f1(a, b, e, x);
 var d1 = hyp2f1(a, b, e + 1.0, x);
 q = a + b + 1.0;
 for (var i = 0; i < aid; i++) {
 r = e - 1.0;
 y = (e * (r - (2.0 * e - q) * x) * d2 + (e - a) * (e - b) * x * d1) / (e * r * s);
 e = r;
 d1 = d2;
 d2 = y;
 }

 return hypdon(y, errObj);
 }

 // negative integer c-a or c-b
 if (negIntCaOrCb) {
 return hypf(a, b, c, d, s, x, errObj);
 }

 return hypok(a, b, c, x, errObj);
};

function hypok(a, b, c, x, errObj) {
 var y = (0, _hyt2f.hyt2f1)(a, b, c, x, errObj);
 return hypdon(y, errObj);
};

function hypdon(y, errObj) {
 if (errObj.val > constants.ETHRESH) {
 // printf( "Estimated err = %.2e\n", err );
 // mtherr("hyp2f1", PLOSS);
 }
 return y;
};

// The transformation for c-a or c-b negative integer
function hypf(a, b, c, d, s, x, errObj) {
 var y = Math.pow(s, d) * (0, _hys2f.hys2f1)(c - a, c - b, c, x, errObj);
 return hypdon(y, errObj);
};

// The alarm exit
function hypdiv() {
 // mtherr("hyp2f1", OVERFLOW);
 return Infinity;
};

exports.hyp2f1 = hyp2f1;
},{"./gamma.js":55,"./hyp2f1/constants.js":58,"./hyp2f1/hys2f1.js":60,"./hyp2f1/hyt2f1.js":61}],58:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
var EPS = exports.EPS = 1.0e-13;
var ETHRESH = exports.ETHRESH = 1.0e-12;
var MAX_ITERATIONS = exports.MAX_ITERATIONS = 10000;
var MACHEP = exports.MACHEP = Number.EPSILON / 2;
},{}],59:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.hyp2f1ra = undefined;

var _constants = require('./constants.js');

var constants = _interopRequireWildcard(_constants);

var _hys2f = require('./hys2f1.js');

function _interopRequireWildcard(obj) { if (obj && obj.__esModule) { return obj; } else { var newObj = {}; if (obj != null) { for (var key in obj) { if (Object.prototype.hasOwnProperty.call(obj, key)) newObj[key] = obj[key]; } } newObj.default = obj; return newObj; } }

function hyp2f1ra(a, b, c, x, lossObj) {
 var da = void 0,
 f2 = void 0,
 f1 = void 0,
 f0 = void 0,
 err = void 0,
 n = void 0;

 // Don't cross c or zero
 if (c < 0 && a <= c || c >= 0 && a >= c) {
 da = Math.round(a - c);
 } else {
 da = Math.round(a);
 }
 var t = a - da;

 lossObj.val = 0;

 // assert(da !== 0);

 if (Math.abs(da) > constants.MAX_ITERATIONS) {
 /* Too expensive to compute this value, so give up */
 // mtherr("hyp2f1", TLOSS);
 lossObj.val = 1.0;
 return NaN;
 }

 if (da < 0) {
 // Recurse down
 f2 = 0;
 f1 = (0, _hys2f.hys2f1)(t, b, c, x, err);
 lossObj.val += err;
 f0 = (0, _hys2f.hys2f1)(t - 1, b, c, x, err);
 lossObj.val += err;
 t -= 1;
 for (n = 1; n < -da; ++n) {
 f2 = f1;
 f1 = f0;
 f0 = -(2 * t - c - t * x + b * x) / (c - t) * f1 - t * (x - 1) / (c - t) * f2;
 t -= 1;
 }
 } else {
 // Recurse up
 f2 = 0;
 f1 = (0, _hys2f.hys2f1)(t, b, c, x, err);
 lossObj.val += err;
 f0 = (0, _hys2f.hys2f1)(t + 1, b, c, x, err);
 lossObj.val += err;
 t += 1;
 for (n = 1; n < da; ++n) {
 f2 = f1;
 f1 = f0;
 f0 = -((2 * t - c - t * x + b * x) * f1 + (c - t) * f2) / (t * (x - 1));
 t += 1;
 }
 }

 return f0;
};

exports.hyp2f1ra = hyp2f1ra;
},{"./constants.js":58,"./hys2f1.js":60}],60:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.hys2f1 = undefined;

var _constants = require('./constants.js');

var constants = _interopRequireWildcard(_constants);

var _hyp2f1ra = require('./hyp2f1ra.js');

function _interopRequireWildcard(obj) { if (obj && obj.__esModule) { return obj; } else { var newObj = {}; if (obj != null) { for (var key in obj) { if (Object.prototype.hasOwnProperty.call(obj, key)) newObj[key] = obj[key]; } } newObj.default = obj; return newObj; } }

function hys2f1(a, b, c, x, lossObj) {
 var intflag = 0;
 var f = void 0;
 // Ensure that |a| > |b| ...
 if (Math.abs(b) > Math.abs(a)) {
 f = b;
 b = a;
 a = f;
 }

 var intB = Math.round(b);

 if (Math.abs(b - intB) < constants.EPS && intB <= 0 && Math.abs(b) < Math.abs(a)) {
 // except when `b` is a smaller negative integer
 f = b;
 b = a;
 a = f;
 intflag = 1;
 }

 if ((Math.abs(a) > Math.abs(c) + 1 || intflag) && Math.abs(c - a) > 2 && Math.abs(a) > 2) {
 // |a| >> |c| implies that large cancellation error is to be expected.
 // We try to reduce it with the recurrence relations
 return (0, _hyp2f1ra.hyp2f1ra)(a, b, c, x, lossObj);
 }

 var i = 0;
 var umax = 0.0;
 f = a;
 var g = b;
 var h = c;
 var s = 1.0;
 var u = 1.0;
 var k = 0.0;
 do {
 if (Math.abs(h) < constants.EPS) {
 lossObj.val = 1.0;
 return Infinity;
 }
 var m = k + 1.0;
 u = u * ((f + k) * (g + k) * x / ((h + k) * m));
 s += u;
 k = Math.abs(u); // remember largest term summed
 if (k > umax) {
 umax = k;
 }
 k = m;
 // should never happen
 if (++i > constants.MAX_ITERATIONS) {
 lossObj.val = 1.0;
 return s;
 }
 } while (s === 0 || Math.abs(u / s) > constants.MACHEP);

 /* return estimated relative error */
 lossObj.val = constants.MACHEP * umax / Math.abs(s) + constants.MACHEP * i;

 return s;
};

exports.hys2f1 = hys2f1;
},{"./constants.js":58,"./hyp2f1ra.js":59}],61:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.hyt2f1 = undefined;

var _slicedToArray = function () { function sliceIterator(arr, i) { var _arr = []; var _n = true; var _d = false; var _e = undefined; try { for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) { _arr.push(_s.value); if (i && _arr.length === i) break; } } catch (err) { _d = true; _e = err; } finally { try { if (!_n && _i["return"]) _i["return"](); } finally { if (_d) throw _e; } } return _arr; } return function (arr, i) { if (Array.isArray(arr)) { return arr; } else if (Symbol.iterator in Object(arr)) { return sliceIterator(arr, i); } else { throw new TypeError("Invalid attempt to destructure non-iterable instance"); } }; }();

var _constants = require('./constants.js');

var constants = _interopRequireWildcard(_constants);

var _gamma = require('../gamma.js');

var _lgam = require('../gamma/lgam.js');

var _psi = require('../psi.js');

var _hys2f = require('./hys2f1.js');

function _interopRequireWildcard(obj) { if (obj && obj.__esModule) { return obj; } else { var newObj = {}; if (obj != null) { for (var key in obj) { if (Object.prototype.hasOwnProperty.call(obj, key)) newObj[key] = obj[key]; } } newObj.default = obj; return newObj; } }

function hyt2f1(a, b, c, x, lossObj) {
 var negIntA = 0;
 var negIntB = 0;
 var e = void 0,
 d1 = void 0,
 d2 = void 0,
 aid = void 0,
 y = void 0,
 t = void 0,
 r = void 0,
 q = void 0,
 sgngam = void 0,
 res = void 0,
 y1 = void 0;

 var aInt = Math.round(a);
 var bInt = Math.round(b);

 if (a <= 0 && Math.abs(a - aInt) < constants.EPS) {
 // a is a negative integer
 negIntA = 1;
 }

 if (b <= 0 && Math.abs(b - bInt) < constants.EPS) {
 // b is a negative integer
 negIntB = 1;
 }

 var errObj = { val: 0.0 };
 var err1Obj = { val: null // note: in original cephes lib err1 is initialized with no value.
 };var s = 1.0 - x;
 if (x < -0.5 && !(negIntA || negIntB)) {
 if (b > a) {
 y = Math.pow(s, -a) * (0, _hys2f.hys2f1)(a, c - b, c, -x / s, errObj);
 } else {
 y = Math.pow(s, -b) * (0, _hys2f.hys2f1)(c - a, b, c, -x / s, errObj);
 }
 return y;
 }

 var d = c - a - b;
 var dInt = Math.round(d);

 if (x > 0.9 && !(negIntA || negIntB)) {
 if (Math.abs(d - dInt) > constants.EPS) {
 // test for integer c-a-b
 // Try the power series first
 y = (0, _hys2f.hys2f1)(a, b, c, x, errObj);
 if (errObj.val < constants.ETHRESH) {
 return y;
 }

 // If power series fails, then apply AMS55 #15.3.6
 q = (0, _hys2f.hys2f1)(a, b, 1.0 - d, s, errObj);
 var sign = 1;
 var w = void 0;
 // lgamSgn(d, &sgngam);

 var _lgamSgn = (0, _lgam.lgamSgn)(d);

 var _lgamSgn2 = _slicedToArray(_lgamSgn, 2);

 w = _lgamSgn2[0];
 sgngam = _lgamSgn2[1];

 sign *= sgngam;

 var _lgamSgn3 = (0, _lgam.lgamSgn)(c - a);

 var _lgamSgn4 = _slicedToArray(_lgamSgn3, 2);

 res = _lgamSgn4[0];
 sgngam = _lgamSgn4[1];

 w -= res;
 sign *= sgngam;

 var _lgamSgn5 = (0, _lgam.lgamSgn)(c - b);

 var _lgamSgn6 = _slicedToArray(_lgamSgn5, 2);

 res = _lgamSgn6[0];
 sgngam = _lgamSgn6[1];

 w -= res;
 sign *= sgngam;
 q *= sign * Math.exp(w);
 r = Math.pow(s, d) * (0, _hys2f.hys2f1)(c - a, c - b, d + 1.0, s, err1Obj);
 sign = 1;

 var _lgamSgn7 = (0, _lgam.lgamSgn)(-d);

 var _lgamSgn8 = _slicedToArray(_lgamSgn7, 2);

 w = _lgamSgn8[0];
 sgngam = _lgamSgn8[1];

 sign *= sgngam;

 var _lgamSgn9 = (0, _lgam.lgamSgn)(a);

 var _lgamSgn10 = _slicedToArray(_lgamSgn9, 2);

 res = _lgamSgn10[0];
 sgngam = _lgamSgn10[1];

 w -= res;
 sign *= sgngam;

 var _lgamSgn11 = (0, _lgam.lgamSgn)(b);

 var _lgamSgn12 = _slicedToArray(_lgamSgn11, 2);

 res = _lgamSgn12[0];
 sgngam = _lgamSgn12[1];

 w -= res;
 sign *= sgngam;
 r *= sign * Math.exp(w);
 y = q + r;

 // estimate cancellation error
 q = Math.abs(q);
 r = Math.abs(r);
 if (q > r) {
 r = q;
 }
 errObj.val += err1Obj.val + constants.MACHEP * r / y;

 y *= (0, _gamma.gamma)(c);
 return y;
 } else {
 /* Psi function expansion, AMS55 #15.3.10, #15.3.11, #15.3.12
 *
 * Although AMS55 does not explicitly state it, this expansion fails
 * for negative integer a or b, since the psi and Gamma functions
 * involved have poles.
 */

 if (dInt >= 0.0) {
 e = d;
 d1 = d;
 d2 = 0.0;
 aid = dInt;
 } else {
 e = -d;
 d1 = 0.0;
 d2 = d;
 aid = -dInt;
 }

 var ax = Math.log(s);

 // sum for t = 0
 y = (0, _psi.psi)(1.0) + (0, _psi.psi)(1.0 + e) - (0, _psi.psi)(a + d1) - (0, _psi.psi)(b + d1) - ax;
 y /= (0, _gamma.gamma)(e + 1.0);

 // Poch for t=1
 var p = (a + d1) * (b + d1) * s / (0, _gamma.gamma)(e + 2.0);
 t = 1.0;
 do {
 r = (0, _psi.psi)(1.0 + t) + (0, _psi.psi)(1.0 + t + e) - (0, _psi.psi)(a + t + d1) - (0, _psi.psi)(b + t + d1) - ax;
 q = p * r;
 y += q;
 p *= s * (a + t + d1) / (t + 1.0);
 p *= (b + t + d1) / (t + 1.0 + e);
 t += 1.0;
 // should never happen
 if (t > constants.MAX_ITERATIONS) {
 // mtherr("hyp2f1", TOOMANY);
 lossObj.val = 1.0;
 return NaN;
 }
 } while (y === 0 || Math.abs(q / y) > constants.EPS);

 if (dInt === 0.0) {
 y *= (0, _gamma.gamma)(c) / ((0, _gamma.gamma)(a) * (0, _gamma.gamma)(b));
 return y;
 }

 y1 = 1.0;

 if (aid === 1) {
 return nosum(a, aid, b, c, d1, d2, dInt, e, q, s, y, y1);
 }

 t = 0.0;
 p = 1.0;
 for (var i = 1; i < aid; i++) {
 r = 1.0 - e + t;
 p *= s * (a + t + d2) * (b + t + d2) / r;
 t += 1.0;
 p /= t;
 y1 += p;
 }
 return nosum(a, aid, b, c, d1, d2, dInt, e, q, s, y, y1);
 }
 }

 // Use defining power series if no special cases
 y = (0, _hys2f.hys2f1)(a, b, c, x, errObj);
 return y;
};

function nosum(a, aid, b, c, d1, d2, dInt, e, q, s, y, y1) {
 var p = (0, _gamma.gamma)(c);
 y1 *= (0, _gamma.gamma)(e) * p / ((0, _gamma.gamma)(a + d1) * (0, _gamma.gamma)(b + d1));

 y *= p / ((0, _gamma.gamma)(a + d2) * (0, _gamma.gamma)(b + d2));
 if ((aid & 1) !== 0) {
 y = -y;
 }
 q = Math.pow(s, dInt);
 if (dInt > 0.0) {
 y *= q;
 } else {
 y1 *= q;
 }

 y += y1;
 return y;
};

exports.hyt2f1 = hyt2f1;
},{"../gamma.js":55,"../gamma/lgam.js":56,"../psi.js":79,"./constants.js":58,"./hys2f1.js":60}],62:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.i0e = exports.i0 = undefined;

var _chbevl = require('./chbevl.js');

/* Chebyshev coefficients for exp(-x) I0(x)
 * in the interval [0,8].
 *
 * lim(x->0){ exp(-x) I0(x) } = 1.
 */
var A = [-4.41534164647933937950E-18, 3.33079451882223809783E-17, -2.43127984654795469359E-16, 1.71539128555513303061E-15, -1.16853328779934516808E-14, 7.67618549860493561688E-14, -4.85644678311192946090E-13, 2.95505266312963983461E-12, -1.72682629144155570723E-11, 9.67580903537323691224E-11, -5.18979560163526290666E-10, 2.65982372468238665035E-9, -1.30002500998624804212E-8, 6.04699502254191894932E-8, -2.67079385394061173391E-7, 1.11738753912010371815E-6, -4.41673835845875056359E-6, 1.64484480707288970893E-5, -5.75419501008210370398E-5, 1.88502885095841655729E-4, -5.76375574538582365885E-4, 1.63947561694133579842E-3, -4.32430999505057594430E-3, 1.05464603945949983183E-2, -2.37374148058994688156E-2, 4.93052842396707084878E-2, -9.49010970480476444210E-2, 1.71620901522208775349E-1, -3.04682672343198398683E-1, 6.76795274409476084995E-1];

/* Chebyshev coefficients for exp(-x) sqrt(x) I0(x)
 * in the inverted interval [8,infinity].
 *
 * lim(x->inf){ exp(-x) sqrt(x) I0(x) } = 1/sqrt(2pi).
 */
/**
 * @file i0.js Modified Bessel function of order zero
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, i0();
 *
 * y = i0( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns modified Bessel function of order zero of the
 * argument.
 *
 * The function is defined as i0(x) = j0( ix ).
 *
 * The range is partitioned into the two intervals [0,8] and
 * (8, infinity). Chebyshev polynomial expansions are employed
 * in each interval.
 *
 *
 *
 * ACCURACY:
 *
 * Relative error:
 * arithmetic domain # trials peak rms
 * IEEE 0,30 30000 5.8e-16 1.4e-16
 *
 *
 * Modified Bessel function of order zero,
 * exponentially scaled
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, i0e();
 *
 * y = i0e( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns exponentially scaled modified Bessel function
 * of order zero of the argument.
 *
 * The function is defined as i0e(x) = exp(-|x|) j0( ix ).
 *
 *
 *
 * ACCURACY:
 *
 * Relative error:
 * arithmetic domain # trials peak rms
 * IEEE 0,30 30000 5.4e-16 1.2e-16
 * See i0().
 *
 * Cephes Math Library Release 2.8: June, 2000
 * Copyright 1984, 1987, 2000 by Stephen L. Moshier
 * Ported to ECMAScript 2018
 * Copyright (c) 2018, Kings Distributed Systems
 *
 * @author KC Erb, kc@kcerb.com
 * @date April 2018
 */
var B = [-7.23318048787475395456E-18, -4.83050448594418207126E-18, 4.46562142029675999901E-17, 3.46122286769746109310E-17, -2.82762398051658348494E-16, -3.42548561967721913462E-16, 1.77256013305652638360E-15, 3.81168066935262242075E-15, -9.55484669882830764870E-15, -4.15056934728722208663E-14, 1.54008621752140982691E-14, 3.85277838274214270114E-13, 7.18012445138366623367E-13, -1.79417853150680611778E-12, -1.32158118404477131188E-11, -3.14991652796324136454E-11, 1.18891471078464383424E-11, 4.94060238822496958910E-10, 3.39623202570838634515E-9, 2.26666899049817806459E-8, 2.04891858946906374183E-7, 2.89137052083475648297E-6, 6.88975834691682398426E-5, 3.36911647825569408990E-3, 8.04490411014108831608E-1];

function i0(x) {
 var y = void 0;
 if (x < 0) x = -x;
 if (x <= 8.0) {
 y = x / 2.0 - 2.0;
 return Math.exp(x) * (0, _chbevl.chbevl)(y, A, 30);
 }
 return Math.exp(x) * (0, _chbevl.chbevl)(32.0 / x - 2.0, B, 25) / Math.sqrt(x);
}

function i0e(x) {
 if (x < 0) x = -x;
 if (x <= 8.0) {
 var y = x / 2.0 - 2.0;
 return (0, _chbevl.chbevl)(y, A, 30);
 }

 return (0, _chbevl.chbevl)(32.0 / x - 2.0, B, 25) / Math.sqrt(x);
}

exports.i0 = i0;
exports.i0e = i0e;
},{"./chbevl.js":50}],63:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.i1e = exports.i1 = undefined;

var _chbevl = require('./chbevl.js');

/* Chebyshev coefficients for exp(-x) I1(x) / x
 * in the interval [0,8].
 *
 * lim(x->0){ exp(-x) I1(x) / x } = 1/2.
 */
var A = [2.77791411276104639959E-18, -2.11142121435816608115E-17, 1.55363195773620046921E-16, -1.10559694773538630805E-15, 7.60068429473540693410E-15, -5.04218550472791168711E-14, 3.22379336594557470981E-13, -1.98397439776494371520E-12, 1.17361862988909016308E-11, -6.66348972350202774223E-11, 3.62559028155211703701E-10, -1.88724975172282928790E-9, 9.38153738649577178388E-9, -4.44505912879632808065E-8, 2.00329475355213526229E-7, -8.56872026469545474066E-7, 3.47025130813767847674E-6, -1.32731636560394358279E-5, 4.78156510755005422638E-5, -1.61760815825896745588E-4, 5.12285956168575772895E-4, -1.51357245063125314899E-3, 4.15642294431288815669E-3, -1.05640848946261981558E-2, 2.47264490306265168283E-2, -5.29459812080949914269E-2, 1.02643658689847095384E-1, -1.76416518357834055153E-1, 2.52587186443633654823E-1];

/* Chebyshev coefficients for exp(-x) sqrt(x) I1(x)
 * in the inverted interval [8,infinity].
 *
 * lim(x->inf){ exp(-x) sqrt(x) I1(x) } = 1/sqrt(2pi).
 */
/**
 * @file i1.js Modified Bessel function of order one
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, i1();
 *
 * y = i1( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns modified Bessel function of order one of the
 * argument.
 *
 * The function is defined as i1(x) = -i j1( ix ).
 *
 * The range is partitioned into the two intervals [0,8] and
 * (8, infinity). Chebyshev polynomial expansions are employed
 * in each interval.
 *
 *
 *
 * ACCURACY:
 *
 * Relative error:
 * arithmetic domain # trials peak rms
 * IEEE 0, 30 30000 1.9e-15 2.1e-16
 *
 *
 * Modified Bessel function of order one,
 * exponentially scaled
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, i1e();
 *
 * y = i1e( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns exponentially scaled modified Bessel function
 * of order one of the argument.
 *
 * The function is defined as i1(x) = -i exp(-|x|) j1( ix ).
 *
 *
 *
 * ACCURACY:
 *
 * Relative error:
 * arithmetic domain # trials peak rms
 * IEEE 0, 30 30000 2.0e-15 2.0e-16
 * See i1().
 *
 * Cephes Math Library Release 2.8: June, 2000
 * Copyright 1985, 1987, 2000 by Stephen L. Moshier
 * Ported to ECMAScript 2018
 * Copyright (c) 2018, Kings Distributed Systems
 *
 * @author KC Erb, kc@kcerb.com
 * @date April 2018
 */
var B = [7.51729631084210481353E-18, 4.41434832307170791151E-18, -4.65030536848935832153E-17, -3.20952592199342395980E-17, 2.96262899764595013876E-16, 3.30820231092092828324E-16, -1.88035477551078244854E-15, -3.81440307243700780478E-15, 1.04202769841288027642E-14, 4.27244001671195135429E-14, -2.10154184277266431302E-14, -4.08355111109219731823E-13, -7.19855177624590851209E-13, 2.03562854414708950722E-12, 1.41258074366137813316E-11, 3.25260358301548823856E-11, -1.89749581235054123450E-11, -5.58974346219658380687E-10, -3.83538038596423702205E-9, -2.63146884688951950684E-8, -2.51223623787020892529E-7, -3.88256480887769039346E-6, -1.10588938762623716291E-4, -9.76109749136146840777E-3, 7.78576235018280120474E-1];

function i1(x) {
 var y = void 0,
 z = void 0;

 z = Math.abs(x);
 if (z <= 8.0) {
 y = z / 2.0 - 2.0;
 z = (0, _chbevl.chbevl)(y, A, 29) * z * Math.exp(z);
 } else {
 z = Math.exp(z) * (0, _chbevl.chbevl)(32.0 / z - 2.0, B, 25) / Math.sqrt(z);
 }
 if (x < 0.0) {
 z = -z;
 }
 return z;
}

function i1e(x) {
 var y = void 0,
 z = void 0;

 z = Math.abs(x);
 if (z <= 8.0) {
 y = z / 2.0 - 2.0;
 z = (0, _chbevl.chbevl)(y, A, 29) * z;
 } else {
 z = (0, _chbevl.chbevl)(32.0 / z - 2.0, B, 25) / Math.sqrt(z);
 }
 if (x < 0.0) {
 z = -z;
 }
 return z;
}

exports.i1 = i1;
exports.i1e = i1e;
},{"./chbevl.js":50}],64:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.iv = undefined;

var _slicedToArray = function () { function sliceIterator(arr, i) { var _arr = []; var _n = true; var _d = false; var _e = undefined; try { for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) { _arr.push(_s.value); if (i && _arr.length === i) break; } } catch (err) { _d = true; _e = err; } finally { try { if (!_n && _i["return"]) _i["return"](); } finally { if (_d) throw _e; } } return _arr; } return function (arr, i) { if (Array.isArray(arr)) { return arr; } else if (Symbol.iterator in Object(arr)) { return sliceIterator(arr, i); } else { throw new TypeError("Invalid attempt to destructure non-iterable instance"); } }; }(); /**
 * @file iv.js Modified Bessel function of noninteger order
 *
 *
 *
 * SYNOPSIS:
 *
 * double v, x, y, iv();
 *
 * y = iv( v, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns modified Bessel function of order v of the
 * argument. If x is negative, v must be integer valued.
 *
 *
 * If x < 0, then v must be an integer.
 * Parts of the code are copyright:
 *
 * Cephes Math Library Release 2.8: June, 2000
 * Copyright 1984, 1987, 1988, 2000 by Stephen L. Moshier
 *
 * And other parts:
 *
 * Copyright (c) 2006 Xiaogang Zhang
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0.
 *
 * Boost Software License - Version 1.0 - August 17th, 2003
 *
 * Permission is hereby granted, free of charge, to any person or
 * organization obtaining a copy of the software and accompanying
 * documentation covered by this license (the "Software") to use, reproduce,
 * display, distribute, execute, and transmit the Software, and to prepare
 * derivative works of the Software, and to permit third-parties to whom the
 * Software is furnished to do so, all subject to the following:
 *
 * The copyright notices in the Software and this entire statement,
 * including the above license grant, this restriction and the following
 * disclaimer, must be included in all copies of the Software, in whole or
 * in part, and all derivative works of the Software, unless such copies or
 * derivative works are solely in the form of machine-executable object code
 * generated by a source language processor.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE AND
 * NON-INFRINGEMENT. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR ANYONE
 * DISTRIBUTING THE SOFTWARE BE LIABLE FOR ANY DAMAGES OR OTHER LIABILITY,
 * WHETHER IN CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 * And the rest are:
 *
 * Copyright (C) 2009 Pauli Virtanen
 * Distributed under the same license as Scipy.
 *
 * And in 2018:
 * Ported to ECMAScript 2018
 * Copyright (c) 2018, Kings Distributed Systems
 *
 * @author KC Erb, kc@kcerb.com
 * @date April 2018
 */


var _ikvTemme3 = require('./iv/ikvTemme.js');

var _ikvAsymptoticUniform3 = require('./iv/ikvAsymptoticUniform.js');

function iv(v, x) {
 var sign = void 0,
 t = void 0,
 ax = void 0,
 res = void 0;

 /* If v is a negative integer, invoke symmetry */
 t = Math.floor(v);
 if (v < 0.0) {
 if (t === v) {
 v = -v; /* symmetry */
 t = -t;
 }
 }
 /* If x is negative, require v to be an integer */
 sign = 1;
 if (x < 0.0) {
 if (t !== v) {
 // mtherr("iv", DOMAIN);
 return NaN;
 }
 if (v !== 2.0 * Math.floor(v / 2.0)) {
 sign = -1;
 }
 }

 /* Avoid logarithm singularity */
 if (x === 0.0) {
 if (v === 0.0) {
 return 1.0;
 }
 if (v < 0.0) {
 // mtherr("iv", OVERFLOW);
 return Infinity;
 } else {
 return 0.0;
 }
 }

 ax = Math.abs(x);
 if (Math.abs(v) > 50) {
 var _ikvAsymptoticUniform = (0, _ikvAsymptoticUniform3.ikvAsymptoticUniform)(v, ax);
 /*
 * Uniform asymptotic expansion for large orders.
 *
 * This appears to overflow slightly later than the Boost
 * implementation of Temme's method.
 */


 var _ikvAsymptoticUniform2 = _slicedToArray(_ikvAsymptoticUniform, 1);

 res = _ikvAsymptoticUniform2[0];
 } else {
 var _ikvTemme = (0, _ikvTemme3.ikvTemme)(v, ax);
 /* Otherwise: Temme's method */


 var _ikvTemme2 = _slicedToArray(_ikvTemme, 1);

 res = _ikvTemme2[0];
 }
 res *= sign;
 return res;
}

exports.iv = iv;
},{"./iv/ikvAsymptoticUniform.js":65,"./iv/ikvTemme.js":66}],65:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.ikvAsymptoticUniform = undefined;

var _constants = require('../constants.js');

var constants = _interopRequireWildcard(_constants);

function _interopRequireWildcard(obj) { if (obj && obj.__esModule) { return obj; } else { var newObj = {}; if (obj != null) { for (var key in obj) { if (Object.prototype.hasOwnProperty.call(obj, key)) newObj[key] = obj[key]; } } newObj.default = obj; return newObj; } }

/*
* Compute Iv, Kv from (AMS5 9.7.7 + 9.7.8), asymptotic expansion for large v
*/
var N_UFACTORS = 11;
var N_UFACTOR_TERMS = 31;

/*
* Uniform asymptotic expansion factors, (AMS5 9.3.9; AMS5 9.3.10)
*/
var asymptoticUfactors = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.20833333333333334, 0.0, 0.125, 0.0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.3342013888888889, 0.0, -0.40104166666666669, 0.0, 0.0703125, 0.0, 0.0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1.0258125964506173, 0.0, 1.8464626736111112, 0.0, -0.89121093750000002, 0.0, 0.0732421875, 0.0, 0.0, 0.0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4.6695844234262474, 0.0, -11.207002616222995, 0.0, 8.78912353515625, 0.0, -2.3640869140624998, 0.0, 0.112152099609375, 0.0, 0.0, 0.0, 0.0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -28.212072558200244, 0.0, 84.636217674600744, 0.0, -91.818241543240035, 0.0, 42.534998745388457, 0.0, -7.3687943594796312, 0.0, 0.22710800170898438, 0.0, 0.0, 0.0, 0.0, 0.0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 212.5701300392171, 0.0, -765.25246814118157, 0.0, 1059.9904525279999, 0.0, -699.57962737613275, 0.0, 218.19051174421159, 0.0, -26.491430486951554, 0.0, 0.57250142097473145, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0, 0, 0, 0, 0, 0, 0, 0, 0, -1919.4576623184068, 0.0, 8061.7221817373083, 0.0, -13586.550006434136, 0.0, 11655.393336864536, 0.0, -5305.6469786134048, 0.0, 1200.9029132163525, 0.0, -108.09091978839464, 0.0, 1.7277275025844574, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0, 0, 0, 0, 0, 0, 20204.291330966149, 0.0, -96980.598388637503, 0.0, 192547.0012325315, 0.0, -203400.17728041555, 0.0, 122200.46498301747, 0.0, -41192.654968897557, 0.0, 7109.5143024893641, 0.0, -493.915304773088, 0.0, 6.074042001273483, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0, 0, 0, -242919.18790055133, 0.0, 1311763.6146629769, 0.0, -2998015.9185381061, 0.0, 3763271.2976564039, 0.0, -2813563.2265865342, 0.0, 1268365.2733216248, 0.0, -331645.17248456361, 0.0, 45218.768981362737, 0.0, -2499.8304818112092, 0.0, 24.380529699556064, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [3284469.8530720375, 0.0, -19706819.11843222, 0.0, 50952602.492664628, 0.0, -74105148.211532637, 0.0, 66344512.274729028, 0.0, -37567176.660763353, 0.0, 13288767.166421819, 0.0, -2785618.1280864552, 0.0, 308186.40461266245, 0.0, -13886.089753717039, 0.0, 110.01714026924674, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]];
function ikvAsymptoticUniform(v, x) {
 // CHECK: I believe the following is equivalent to what we had before ....
 var iValue = 0;
 var kValue = null;
 var iPrefactor = void 0,
 kPrefactor = void 0,
 t = void 0,
 t2 = void 0,
 eta = void 0,
 z = void 0,
 iSum = void 0,
 kSUm = void 0,
 term = void 0,
 divisor = void 0,
 k = void 0,
 n = void 0;
 var sign = 1;

 if (v < 0) {
 /* Negative v; compute I_{-v} and K_{-v} and use (AMS 9.6.2) */
 sign = -1;
 v = -v;
 }

 z = x / v;
 t = 1 / Math.sqrt(1 + z * z);
 t2 = t * t;
 eta = Math.sqrt(1 + z * z) + Math.log(z / (1 + 1 / t));

 iPrefactor = Math.sqrt(t / (2 * Math.PI * v)) * Math.exp(v * eta);
 iSum = 1.0;

 kPrefactor = Math.sqrt(Math.PI * t / (2 * v)) * Math.exp(-v * eta);
 kSUm = 1.0;

 divisor = v;
 for (n = 1; n < N_UFACTORS; ++n) {
 /*
 * Evaluate u_k(t) with Horner's scheme;
 * (using the knowledge about which coefficients are zero)
 */
 term = 0;
 for (k = N_UFACTOR_TERMS - 1 - 3 * n; k < N_UFACTOR_TERMS - n; k += 2) {
 term *= t2;
 term += asymptoticUfactors[n][k];
 }
 for (k = 1; k < n; k += 2) {
 term *= t2;
 }
 if (n % 2 === 1) {
 term *= t;
 }

 /* Sum terms */
 term /= divisor;
 iSum += term;
 kSUm += n % 2 === 0 ? term : -term;

 /* Check convergence */
 if (Math.abs(term) < constants.MACHEP) {
 break;
 }

 divisor *= v;
 }

 if (Math.abs(term) > 1e-3 * Math.abs(iSum)) {
 /* Didn't converge */
 // mtherr("ikvAsymptoticUniform", TLOSS);
 }
 if (Math.abs(term) > constants.MACHEP * Math.abs(iSum)) {
 /* Some precision lost */
 // mtherr("ikvAsymptoticUniform", PLOSS);
 }

 if (kValue !== null) {
 /* symmetric in v */
 kValue = kPrefactor * kSUm;
 }

 if (iValue !== null) {
 if (sign === 1) {
 iValue = iPrefactor * iSum;
 } else {
 /* (AMS 9.6.2) */
 iValue = iPrefactor * iSum + 2 / Math.PI * Math.sin(Math.PI * v) * kPrefactor * kSUm;
 }
 }

 return [iValue, kValue];
}

exports.ikvAsymptoticUniform = ikvAsymptoticUniform;
},{"../constants.js":51}],66:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.ikvTemme = undefined;

var _slicedToArray = function () { function sliceIterator(arr, i) { var _arr = []; var _n = true; var _d = false; var _e = undefined; try { for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) { _arr.push(_s.value); if (i && _arr.length === i) break; } } catch (err) { _d = true; _e = err; } finally { try { if (!_n && _i["return"]) _i["return"](); } finally { if (_d) throw _e; } } return _arr; } return function (arr, i) { if (Array.isArray(arr)) { return arr; } else if (Symbol.iterator in Object(arr)) { return sliceIterator(arr, i); } else { throw new TypeError("Invalid attempt to destructure non-iterable instance"); } }; }();

var _constants = require('../constants.js');

var constants = _interopRequireWildcard(_constants);

var _gamma = require('../gamma.js');

function _interopRequireWildcard(obj) { if (obj && obj.__esModule) { return obj; } else { var newObj = {}; if (obj != null) { for (var key in obj) { if (Object.prototype.hasOwnProperty.call(obj, key)) newObj[key] = obj[key]; } } newObj.default = obj; return newObj; } }

/*
* Compute I(v, x) and K(v, x) simultaneously by Temme's method, see
* Temme, Journal of Computational Physics, vol 19, 324 (1975)
*/
var NEED_I = 0x1;
var NEED_K = 0x2;

function ikvTemme(v, x) {
 var ivP = 0;
 var kvP = null;
 var u = void 0,
 iV = void 0,
 kv = void 0,
 kv1 = void 0,
 ku = void 0,
 ku1 = void 0,
 fv = void 0,
 capitalW = void 0,
 current = void 0,
 prev = void 0,
 next = void 0,
 n = void 0,
 k = void 0,
 kind = void 0,
 z = void 0,
 lim = void 0;
 var reflect = 0;
 kind = 0;
 if (ivP !== null) {
 kind = kind | NEED_I;
 }
 if (kvP !== null) {
 kind = kind | NEED_K;
 }

 if (v < 0) {
 reflect = 1;
 v = -v; /* v is non-negative from here */
 kind = kind | NEED_K;
 }
 n = Math.round(v);
 u = v - n; /* -1/2 <= u < 1/2 */

 if (x < 0) {
 if (ivP !== null) {
 ivP = NaN;
 }
 if (kvP !== null) {
 kvP = NaN;
 }
 // mtherr("ikv_temme", DOMAIN);
 return [ivP, kvP];
 }
 if (x === 0) {
 iV = v === 0 ? 1 : 0;
 if (kind & NEED_K) {
 // mtherr("ikv_temme", OVERFLOW);
 kv = Infinity;
 } else {
 kv = NaN; /* any value will do */
 }

 if (reflect && kind & NEED_I) {
 z = u + n % 2;

 iV = Math.sin(Math.PI * z) === 0 ? iV : Infinity;
 if (iV === Infinity || iV === -Infinity) {
 // mtherr("ikv_temme", OVERFLOW);
 }
 }

 if (ivP !== null) {
 ivP = iV;
 }
 if (kvP !== null) {
 kvP = kv;
 }
 return [ivP, kvP];
 }

 /* x is positive until reflection */
 capitalW = 1 / x; /* Wronskian */
 if (x <= 2) {
 /* Temme series */
 var _temmeIkSeries = temmeIkSeries(u, x); /* x in (0, 2] */


 var _temmeIkSeries2 = _slicedToArray(_temmeIkSeries, 2);

 ku = _temmeIkSeries2[0];
 ku1 = _temmeIkSeries2[1];
 } else {
 /* continued fraction cF2Ik */
 var _cF2Ik = cF2Ik(u, x, ku, ku1); /* x in (2, \infty) */


 var _cF2Ik2 = _slicedToArray(_cF2Ik, 2);

 ku = _cF2Ik2[0];
 ku1 = _cF2Ik2[1];
 }
 prev = ku;
 current = ku1;
 for (k = 1; k <= n; k++) {
 /* forward recurrence for K */
 next = 2 * (u + k) * current / x + prev;
 prev = current;
 current = next;
 }
 kv = prev;
 kv1 = current;
 if (kind & NEED_I) {
 lim = (4 * v * v + 10) / (8 * x);
 lim *= lim;
 lim *= lim;
 lim /= 24;
 if (lim < constants.MACHEP * 10 && x > 100) {
 /*
 * x is huge compared to v, CF1 may be very slow
 * to converge so use asymptotic expansion for large
 * x case instead. Note that the asymptotic expansion
 * isn't very accurate - so it's deliberately very hard
 * to get here - probably we're going to overflow:
 */
 iV = ivAsymptotic(v, x);
 } else {
 fv = cF1Ik(v, x); /* continued fraction cF1Ik */
 iV = capitalW / (kv * fv + kv1); /* Wronskian relation */
 }
 } else {
 iV = NaN; /* any value will do */
 }

 if (reflect) {
 z = u + n % 2;

 if (ivP !== null) {
 ivP = iV + 2 / Math.PI * Math.sin(Math.PI * z) * kv; /* reflection formula */
 }
 if (kvP !== null) {
 kvP = kv;
 }
 } else {
 if (ivP !== null) {
 ivP = iV;
 }
 if (kvP !== null) {
 kvP = kv;
 }
 }
 return [ivP, kvP];
}

/*
* Compute iV from (AMS5 9.7.1), asymptotic expansion for large |z|
* iV ~ exp(x)/Math.sqrt(2 pi x) ( 1 + (4*v*v-1)/8x + (4*v*v-1)(4*v*v-9)/8x/2! + ...)
*/
function ivAsymptotic(v, x) {
 var mu = void 0,
 sum = void 0,
 term = void 0,
 prefactor = void 0,
 factor = void 0,
 k = void 0;
 prefactor = Math.exp(x) / Math.sqrt(2 * Math.PI * x);

 if (prefactor === Infinity) {
 return prefactor;
 }

 mu = 4 * v * v;
 sum = 1.0;
 term = 1.0;
 k = 1;

 do {
 factor = (mu - (2 * k - 1) * (2 * k - 1)) / (8 * x) / k;
 if (k > 100) {
 /* didn't converge */
 // mtherr("iv(ivAsymptotic)", TLOSS);
 break;
 }
 term *= -factor;
 sum += term;
 ++k;
 } while (Math.abs(term) > constants.MACHEP * Math.abs(sum));
 return sum * prefactor;
}

/*
* Modified Bessel functions of the first and second kind of fractional order
*
* Calculate K(v, x) and K(v+1, x) by method analogous to
* Temme, Journal of Computational Physics, vol 21, 343 (1976)
*/
function temmeIkSeries(v, x) {
 var f = void 0,
 h = void 0,
 p = void 0,
 q = void 0,
 coef = void 0,
 sum = void 0,
 sum1 = void 0,
 tolerance = void 0,
 a = void 0,
 b = void 0,
 c = void 0,
 d = void 0,
 sigma = void 0,
 gamma1 = void 0,
 gamma2 = void 0,
 k = void 0,
 gp = void 0,
 gm = void 0,
 kres = void 0,
 kres1 = void 0;

 /*
 * |x| <= 2, Temme series converge rapidly
 * |x| > 2, the larger the |x|, the slower the convergence
 */
 boostAssert(Math.abs(x) <= 2, '|x| > 2, the larger the |x|, the slower the convergence');
 boostAssert(Math.abs(v) <= 0.5, '|v| > 0.5');

 gp = (0, _gamma.gamma)(v + 1) - 1;
 gm = (0, _gamma.gamma)(-v + 1) - 1;

 a = Math.log(x / 2);
 b = Math.exp(v * a);
 sigma = -a * v;
 c = Math.abs(v) < constants.MACHEP ? 1 : Math.sin(Math.PI * v) / (v * Math.PI);
 d = Math.abs(sigma) < constants.MACHEP ? 1 : Math.sinh(sigma) / sigma;
 gamma1 = Math.abs(v) < constants.MACHEP ? -constants.EULER : 0.5 / v * (gp - gm) * c;
 gamma2 = (2 + gp + gm) * c / 2;

 /* initial values */
 p = (gp + 1) / (2 * b);
 q = (1 + gm) * b / 2;
 f = (Math.cosh(sigma) * gamma1 + d * -a * gamma2) / c;
 h = p;
 coef = 1;
 sum = coef * f;
 sum1 = coef * h;

 /* series summation */
 tolerance = constants.MACHEP;
 for (k = 1; k < constants.MAXITER; k++) {
 f = (k * f + p + q) / (k * k - v * v);
 p /= k - v;
 q /= k + v;
 h = p - k * f;
 coef *= x * x / (4 * k);
 sum += coef * f;
 sum1 += coef * h;
 if (Math.abs(coef * f) < Math.abs(sum) * tolerance) {
 break;
 }
 }
 if (k === constants.MAXITER) {
 // mtherr("ikvTemme(temmeIkSeries)", TLOSS);
 }

 kres = sum;
 kres1 = 2 * sum1 / x;

 return [kres, kres1];
}

/* Evaluate continued fraction fv = I_(v+1) / I_v, derived from
* Abramowitz and Stegun, Handbook of Mathematical Functions, 1972, 9.1.73 */
function cF1Ik(v, x) {
 var capitalC = void 0,
 capitalD = void 0,
 f = void 0,
 a = void 0,
 b = void 0,
 delta = void 0,
 tiny = void 0,
 tolerance = void 0,
 fv = void 0,
 k = void 0;
 /*
 * |x| <= |v|, cF1Ik converges rapidly
 * |x| > |v|, cF1Ik needs O(|x|) iterations to converge
 */

 /*
 * modified Lentz's method, see
 * Lentz, Applied Optics, vol 15, 668 (1976)
 */
 tolerance = 2 * constants.MACHEP;
 tiny = 1 / Math.sqrt(constants.DBL_MAX);
 capitalC = f = tiny; /* b0 = 0, replace with tiny */
 capitalD = 0;
 for (k = 1; k < constants.MAXITER; k++) {
 a = 1;
 b = 2 * (v + k) / x;
 capitalC = b + a / capitalC;
 capitalD = b + a * capitalD;
 if (capitalC === 0) {
 capitalC = tiny;
 }
 if (capitalD === 0) {
 capitalD = tiny;
 }
 capitalD = 1 / capitalD;
 delta = capitalC * capitalD;
 f *= delta;
 if (Math.abs(delta - 1) <= tolerance) {
 break;
 }
 }
 if (k === constants.MAXITER) {
 // mtherr("ikvTemme(cF1Ik)", TLOSS);
 }

 fv = f;

 return fv;
}

/*
* Calculate K(v, x) and K(v+1, x) by evaluating continued fraction
* z1 / z0 = U(v+1.5, 2v+1, 2x) / U(v+0.5, 2v+1, 2x), see
* Thompson and Barnett, Computer Physics Communications, vol 47, 245 (1987)
*/
function cF2Ik(v, x) {
 var capitalS = void 0,
 capitalC = void 0,
 capitalQ = void 0,
 capitalD = void 0,
 f = void 0,
 a = void 0,
 b = void 0,
 q = void 0,
 delta = void 0,
 tolerance = void 0,
 current = void 0,
 prev = void 0,
 kvRes = void 0,
 kv1Res = void 0,
 k = void 0;

 /*
 * |x| >= |v|, cF2Ik converges rapidly
 * |x| -> 0, cF2Ik fails to converge
 */
 boostAssert(Math.abs(x) > 1, 'cF2Ik fails to converge since: |x| <= 1');

 /*
 * Steed's algorithm, see Thompson and Barnett,
 * Journal of Computational Physics, vol 64, 490 (1986)
 */
 tolerance = constants.MACHEP;
 a = v * v - 0.25;
 b = 2 * (x + 1); /* b1 */
 capitalD = 1 / b; /* D1 = 1 / b1 */
 f = delta = capitalD; /* f1 = delta1 = D1, coincidence */
 prev = 0; /* q0 */
 current = 1; /* q1 */
 capitalQ = capitalC = -a; /* Q1 = C1 because q1 = 1 */
 capitalS = 1 + capitalQ * delta; /* S1 */
 for (k = 2; k < constants.MAXITER; k++) {
 /* starting from 2 */
 /* continued fraction f = z1 / z0 */
 a -= 2 * (k - 1);
 b += 2;
 capitalD = 1 / (b + a * capitalD);
 delta *= b * capitalD - 1;
 f += delta;

 /* series summation capitalS = 1 + \sum_{n=1}^{\infty} C_n * z_n / z_0 */
 q = (prev - (b - 2) * current) / a;
 prev = current;
 current = q; /* forward recurrence for q */
 capitalC *= -a / k;
 capitalQ += capitalC * q;
 capitalS += capitalQ * delta;

 /* capitalS converges slower than f */
 if (Math.abs(capitalQ * delta) < Math.abs(capitalS) * tolerance) {
 break;
 }
 }
 if (k === constants.MAXITER) {
 // mtherr("ikvTemme(cF2Ik)", TLOSS);
 }

 kvRes = Math.sqrt(Math.PI / (2 * x)) * Math.exp(-x) / capitalS;
 kv1Res = kvRes * (0.5 + v + x + (v * v - 0.25) * f) / x;

 return [kvRes, kv1Res];
}

function boostAssert(expr, msg) {
 if (expr) {
 // do nothing
 } else {
 throw new Error(msg);
 }
}
exports.ikvTemme = ikvTemme;
},{"../constants.js":51,"../gamma.js":55}],67:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.y0 = exports.j0 = undefined;

var _constants = require('./constants.js');

var constants = _interopRequireWildcard(_constants);

var _polevl = require('./polevl.js');

function _interopRequireWildcard(obj) { if (obj && obj.__esModule) { return obj; } else { var newObj = {}; if (obj != null) { for (var key in obj) { if (Object.prototype.hasOwnProperty.call(obj, key)) newObj[key] = obj[key]; } } newObj.default = obj; return newObj; } }

/**
 * @file j0.js Bessel function of order zero
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, j0();
 *
 * y = j0( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of order zero of the argument.
 *
 * The domain is divided into the intervals [0, 5] and
 * (5, infinity). In the first interval the following rational
 * approximation is used:
 *
 *
 * 2 2
 * (w - r ) (w - r ) P (w) / Q (w)
 * 1 2 3 8
 *
 * 2
 * where w = x and the two r's are zeros of the function.
 *
 * In the second interval, the Hankel asymptotic expansion
 * is employed with two rational functions of degree 6/6
 * and 7/7.
 *
 *
 *
 * ACCURACY:
 *
 * Absolute error:
 * arithmetic domain # trials peak rms
 * IEEE 0, 30 60000 4.2e-16 1.1e-16
 *
 *
 * Bessel function of the second kind, order zero
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, y0();
 *
 * y = y0( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of the second kind, of order
 * zero, of the argument.
 *
 * The domain is divided into the intervals [0, 5] and
 * (5, infinity). In the first interval a rational approximation
 * R(x) is employed to compute
 * y0(x) = R(x) + 2 * log(x) * j0(x) / NPY_PI.
 * Thus a call to j0() is required.
 *
 * In the second interval, the Hankel asymptotic expansion
 * is employed with two rational functions of degree 6/6
 * and 7/7.
 *
 *
 *
 * ACCURACY:
 *
 * Absolute error, when y0(x) < 1; else relative error:
 *
 * arithmetic domain # trials peak rms
 * IEEE 0, 30 30000 1.3e-15 1.6e-16
 *
 * Cephes Math Library Release 2.8: June, 2000
 * Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
 * Ported to ECMAScript 2018
 * Copyright (c) 2018, Kings Distributed Systems
 *
 * @author KC Erb, kc@kcerb.com
 * @date April 2018
 */

/* Note: all coefficients satisfy the relative error criterion
 * except YP, YQ which are designed for absolute error. */
var PP = [7.96936729297347051624E-4, 8.28352392107440799803E-2, 1.23953371646414299388E0, 5.44725003058768775090E0, 8.74716500199817011941E0, 5.30324038235394892183E0, 9.99999999999999997821E-1];

var PQ = [9.24408810558863637013E-4, 8.56288474354474431428E-2, 1.25352743901058953537E0, 5.47097740330417105182E0, 8.76190883237069594232E0, 5.30605288235394617618E0, 1.00000000000000000218E0];

var QP = [-1.13663838898469149931E-2, -1.28252718670509318512E0, -1.95539544257735972385E1, -9.32060152123768231369E1, -1.77681167980488050595E2, -1.47077505154951170175E2, -5.14105326766599330220E1, -6.05014350600728481186E0];

var QQ = [
/* 1.00000000000000000000E0, */
6.43178256118178023184E1, 8.56430025976980587198E2, 3.88240183605401609683E3, 7.24046774195652478189E3, 5.93072701187316984827E3, 2.06209331660327847417E3, 2.42005740240291393179E2];

var YP = [1.55924367855235737965E4, -1.46639295903971606143E7, 5.43526477051876500413E9, -9.82136065717911466409E11, 8.75906394395366999549E13, -3.46628303384729719441E15, 4.42733268572569800351E16, -1.84950800436986690637E16];

var YQ = [
/* 1.00000000000000000000E0, */
1.04128353664259848412E3, 6.26107330137134956842E5, 2.68919633393814121987E8, 8.64002487103935000337E10, 2.02979612750105546709E13, 3.17157752842975028269E15, 2.50596256172653059228E17];

/* 5.783185962946784521175995758455807035071 */
var DR1 = 5.78318596294678452118E0;

/* 30.47126234366208639907816317502275584842 */
var DR2 = 3.04712623436620863991E1;

var RP = [-4.79443220978201773821E9, 1.95617491946556577543E12, -2.49248344360967716204E14, 9.70862251047306323952E15];

var RQ = [
/* 1.00000000000000000000E0, */
4.99563147152651017219E2, 1.73785401676374683123E5, 4.84409658339962045305E7, 1.11855537045356834862E10, 2.11277520115489217587E12, 3.10518229857422583814E14, 3.18121955943204943306E16, 1.71086294081043136091E18];

function j0(x) {
 var w = void 0,
 z = void 0,
 p = void 0,
 q = void 0,
 xn = void 0;

 if (x < 0) x = -x;

 if (x <= 5.0) {
 z = x * x;
 if (x < 1.0e-5) return 1.0 - z / 4.0;

 p = (z - DR1) * (z - DR2);
 p = p * (0, _polevl.polevl)(z, RP, 3) / (0, _polevl.p1evl)(z, RQ, 8);
 return p;
 }

 w = 5.0 / x;
 q = 25.0 / (x * x);
 p = (0, _polevl.polevl)(q, PP, 6) / (0, _polevl.polevl)(q, PQ, 6);
 q = (0, _polevl.polevl)(q, QP, 7) / (0, _polevl.p1evl)(q, QQ, 7);
 xn = x - constants.NPY_PI_4;
 p = p * Math.cos(xn) - w * q * Math.sin(xn);
 return p * constants.SQ2OPI / Math.sqrt(x);
}

/* y0() 2 */
/* Bessel function of second kind, order zero */

/* Rational approximation coefficients YP[], YQ[] are used here.
 * The function computed is y0(x) - 2 * log(x) * j0(x) / NPY_PI,
 * whose value at x = 0 is 2 * ( log(0.5) + EUL ) / NPY_PI
 * = 0.073804295108687225.
 */

function y0(x) {
 var w = void 0,
 z = void 0,
 p = void 0,
 q = void 0,
 xn = void 0;

 if (x <= 5.0) {
 if (x === 0.0) {
 // mtherr('y0', SING);
 return -Infinity;
 } else if (x < 0.0) {
 // mtherr('y0', DOMAIN);
 return NaN;
 }

 z = x * x;
 w = (0, _polevl.polevl)(z, YP, 7) / (0, _polevl.p1evl)(z, YQ, 7);
 w += constants.NPY_2_PI * Math.log(x) * j0(x);
 return w;
 }

 w = 5.0 / x;
 z = 25.0 / (x * x);
 p = (0, _polevl.polevl)(z, PP, 6) / (0, _polevl.polevl)(z, PQ, 6);
 q = (0, _polevl.polevl)(z, QP, 7) / (0, _polevl.p1evl)(z, QQ, 7);
 xn = x - constants.NPY_PI_4;
 p = p * Math.sin(xn) + w * q * Math.cos(xn);
 return p * constants.SQ2OPI / Math.sqrt(x);
}

exports.j0 = j0;
exports.y0 = y0;
},{"./constants.js":51,"./polevl.js":78}],68:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.y1 = exports.j1 = undefined;

var _constants = require('./constants.js');

var constants = _interopRequireWildcard(_constants);

var _polevl = require('./polevl.js');

function _interopRequireWildcard(obj) { if (obj && obj.__esModule) { return obj; } else { var newObj = {}; if (obj != null) { for (var key in obj) { if (Object.prototype.hasOwnProperty.call(obj, key)) newObj[key] = obj[key]; } } newObj.default = obj; return newObj; } }

/**
 * @file j1.js Bessel function of order one
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, j1();
 *
 * y = j1( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of order one of the argument.
 *
 * The domain is divided into the intervals [0, 8] and
 * (8, infinity). In the first interval a 24 term Chebyshev
 * expansion is used. In the second, the asymptotic
 * trigonometric representation is employed using two
 * rational functions of degree 5/5.
 *
 *
 *
 * ACCURACY:
 *
 * Absolute error:
 * arithmetic domain # trials peak rms
 * IEEE 0, 30 30000 2.6e-16 1.1e-16
 *
 *
 * Bessel function of second kind of order one
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, y1();
 *
 * y = y1( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of the second kind of order one
 * of the argument.
 *
 * The domain is divided into the intervals [0, 8] and
 * (8, infinity). In the first interval a 25 term Chebyshev
 * expansion is used, and a call to j1() is required.
 * In the second, the asymptotic trigonometric representation
 * is employed using two rational functions of degree 5/5.
 *
 *
 *
 * ACCURACY:
 *
 * Absolute error:
 * arithmetic domain # trials peak rms
 * IEEE 0, 30 30000 1.0e-15 1.3e-16
 *
 * (error criterion relative when |y1| > 1).
 *
 * Cephes Math Library Release 2.8: June, 2000
 * Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
 * Ported to ECMAScript 2018
 * Copyright (c) 2018, Kings Distributed Systems
 *
 * @author KC Erb, kc@kcerb.com
 * @date April 2018
 */
var RP = [-8.99971225705559398224E8, 4.52228297998194034323E11, -7.27494245221818276015E13, 3.68295732863852883286E15];

var RQ = [
/* 1.00000000000000000000E0, */
6.20836478118054335476E2, 2.56987256757748830383E5, 8.35146791431949253037E7, 2.21511595479792499675E10, 4.74914122079991414898E12, 7.84369607876235854894E14, 8.95222336184627338078E16, 5.32278620332680085395E18];

var PP = [7.62125616208173112003E-4, 7.31397056940917570436E-2, 1.12719608129684925192E0, 5.11207951146807644818E0, 8.42404590141772420927E0, 5.21451598682361504063E0, 1.00000000000000000254E0];

var PQ = [5.71323128072548699714E-4, 6.88455908754495404082E-2, 1.10514232634061696926E0, 5.07386386128601488557E0, 8.39985554327604159757E0, 5.20982848682361821619E0, 9.99999999999999997461E-1];

var QP = [5.10862594750176621635E-2, 4.98213872951233449420E0, 7.58238284132545283818E1, 3.66779609360150777800E2, 7.10856304998926107277E2, 5.97489612400613639965E2, 2.11688757100572135698E2, 2.52070205858023719784E1];

var QQ = [
/* 1.00000000000000000000E0, */
7.42373277035675149943E1, 1.05644886038262816351E3, 4.98641058337653607651E3, 9.56231892404756170795E3, 7.99704160447350683650E3, 2.82619278517639096600E3, 3.36093607810698293419E2];

var YP = [1.26320474790178026440E9, -6.47355876379160291031E11, 1.14509511541823727583E14, -8.12770255501325109621E15, 2.02439475713594898196E17, -7.78877196265950026825E17];

var YQ = [
/* 1.00000000000000000000E0, */
5.94301592346128195359E2, 2.35564092943068577943E5, 7.34811944459721705660E7, 1.87601316108706159478E10, 3.88231277496238566008E12, 6.20557727146953693363E14, 6.87141087355300489866E16, 3.97270608116560655612E18];

var Z1 = 1.46819706421238932572E1;
var Z2 = 4.92184563216946036703E1;

function j1(x) {
 var w = void 0,
 z = void 0,
 p = void 0,
 q = void 0,
 xn = void 0;

 w = x;
 if (x < 0) return -j1(-x);

 if (w <= 5.0) {
 z = x * x;
 w = (0, _polevl.polevl)(z, RP, 3) / (0, _polevl.p1evl)(z, RQ, 8);
 w = w * x * (z - Z1) * (z - Z2);
 return w;
 }

 w = 5.0 / x;
 z = w * w;
 p = (0, _polevl.polevl)(z, PP, 6) / (0, _polevl.polevl)(z, PQ, 6);
 q = (0, _polevl.polevl)(z, QP, 7) / (0, _polevl.p1evl)(z, QQ, 7);
 xn = x - constants.THPIO4;
 p = p * Math.cos(xn) - w * q * Math.sin(xn);
 return p * constants.SQ2OPI / Math.sqrt(x);
}

function y1(x) {
 var w = void 0,
 z = void 0,
 p = void 0,
 q = void 0,
 xn = void 0;

 if (x <= 5.0) {
 if (x === 0.0) {
 // mtherr("y1", SING);
 return -Infinity;
 } else if (x <= 0.0) {
 // mtherr("y1", DOMAIN);
 return NaN;
 }
 z = x * x;
 w = x * ((0, _polevl.polevl)(z, YP, 5) / (0, _polevl.p1evl)(z, YQ, 8));
 w += constants.NPY_2_PI * (j1(x) * Math.log(x) - 1.0 / x);
 return w;
 }

 w = 5.0 / x;
 z = w * w;
 p = (0, _polevl.polevl)(z, PP, 6) / (0, _polevl.polevl)(z, PQ, 6);
 q = (0, _polevl.polevl)(z, QP, 7) / (0, _polevl.p1evl)(z, QQ, 7);
 xn = x - constants.THPIO4;
 p = p * Math.sin(xn) + w * q * Math.cos(xn);
 return p * constants.SQ2OPI / Math.sqrt(x);
}

exports.j1 = j1;
exports.y1 = y1;
},{"./constants.js":51,"./polevl.js":78}],69:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.jv = undefined;

var _slicedToArray = function () { function sliceIterator(arr, i) { var _arr = []; var _n = true; var _d = false; var _e = undefined; try { for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) { _arr.push(_s.value); if (i && _arr.length === i) break; } } catch (err) { _d = true; _e = err; } finally { try { if (!_n && _i["return"]) _i["return"](); } finally { if (_d) throw _e; } } return _arr; } return function (arr, i) { if (Array.isArray(arr)) { return arr; } else if (Symbol.iterator in Object(arr)) { return sliceIterator(arr, i); } else { throw new TypeError("Invalid attempt to destructure non-iterable instance"); } }; }(); /**
 * @file jv.js Bessel function of noninteger order
 *
 *
 *
 * SYNOPSIS:
 *
 * double v, x, y, jv();
 *
 * y = jv( v, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of order v of the argument,
 * where v is real. Negative x is allowed if v is an integer.
 *
 * Several expansions are included: the ascending power
 * series, the Hankel expansion, and two transitional
 * expansions for large v. If v is not too large, it
 * is reduced by recurrence to a region of best accuracy.
 * The transitional expansions give 12D accuracy for v > 500.
 *
 *
 *
 * ACCURACY:
 * Results for integer v are indicated by *, where x and v
 * both vary from -125 to +125. Otherwise,
 * x ranges from 0 to 125, v ranges as indicated by "domain."
 * Error criterion is absolute, except relative when |jv()| > 1.
 *
 * arithmetic v domain x domain # trials peak rms
 * IEEE 0,125 0,125 100000 4.6e-15 2.2e-16
 * IEEE -125,0 0,125 40000 5.4e-11 3.7e-13
 * IEEE 0,500 0,500 20000 4.4e-15 4.0e-16
 * Integer v:
 * IEEE -125,125 -125,125 50000 3.5e-15* 1.9e-16*
 *
 * Cephes Math Library Release 2.8: June, 2000
 * Copyright 1984, 1987, 1989, 1992, 2000 by Stephen L. Moshier
 * Ported to ECMAScript 2018
 * Copyright (c) 2018, Kings Distributed Systems
 *
 * @author KC Erb, kc@kcerb.com
 * @date April 2018
 */


var _constants = require('./constants.js');

var constants = _interopRequireWildcard(_constants);

var _j = require('./j0.js');

var _j2 = require('./j1.js');

var _gamma = require('./gamma.js');

var _jvs = require('./jv/jvs.js');

var _hankel = require('./jv/hankel.js');

var _recur9 = require('./jv/recur.js');

var _jnx = require('./jv/jnx.js');

function _interopRequireWildcard(obj) { if (obj && obj.__esModule) { return obj; } else { var newObj = {}; if (obj != null) { for (var key in obj) { if (Object.prototype.hasOwnProperty.call(obj, key)) newObj[key] = obj[key]; } } newObj.default = obj; return newObj; } }

function jv(n, x) {
 var k = void 0,
 q = void 0,
 t = void 0,
 y = void 0,
 an = void 0,
 i = void 0,
 sign = void 0,
 nint = void 0,
 condition = void 0;

 nint = 0; /* Flag for integer n */
 sign = 1; /* Flag for sign inversion */
 an = Math.abs(n);
 y = Math.floor(an);
 if (y === an) {
 nint = 1;
 i = an - 16384.0 * Math.floor(an / 16384.0);
 if (n < 0.0) {
 if (i & 1) {
 sign = -sign;
 }
 n = an;
 }
 if (x < 0.0) {
 if (i & 1) {
 sign = -sign;
 }
 x = -x;
 }
 if (n === 0.0) {
 return (0, _j.j0)(x);
 }
 if (n === 1.0) {
 return sign * (0, _j2.j1)(x);
 }
 }

 if (x < 0.0 && y !== an) {
 // mtherr("Jv", DOMAIN);
 y = NaN;
 return sign * y;
 }

 if (x === 0 && n < 0 && !nint) {
 // mtherr("Jv", OVERFLOW);
 return Infinity / (0, _gamma.gamma)(n + 1);
 }

 y = Math.abs(x);

 if (y * y < Math.abs(n + 1) * constants.MACHEP) {
 return Math.pow(0.5 * x, n) / (0, _gamma.gamma)(n + 1);
 }

 k = 3.6 * Math.sqrt(y);
 t = 3.6 * Math.sqrt(an);

 if (y < t && an > 21.0) return sign * (0, _jvs.jvs)(n, x);
 if (an < k && y > 21.0) return sign * (0, _hankel.hankel)(n, x);

 if (an < 500.0) {
 /* Note: if x is too large, the continued fraction will fail; but then the
 * Hankel expansion can be used. */
 if (nint !== 0) {
 k = 0.0;

 var _recur = (0, _recur9.recur)(n, x, k, 1);

 var _recur2 = _slicedToArray(_recur, 3);

 q = _recur2[0];
 n = _recur2[1];
 k = _recur2[2];

 if (k === 0.0) {
 y = (0, _j.j0)(x) / q;
 return sign * y;
 }
 if (k === 1.0) {
 y = (0, _j2.j1)(x) / q;
 return sign * y;
 }
 }

 condition = an > 2.0 * y || n >= 0.0 && n < 20.0 && y > 6.0 && y < 20.0;

 if (condition) {
 /* Recur backwards from a larger value of n */
 k = n;

 y = y + an + 1.0;
 if (y < 30.0) y = 30.0;
 y = n + Math.floor(y - n);

 var _recur3 = (0, _recur9.recur)(y, x, k, 0);

 var _recur4 = _slicedToArray(_recur3, 3);

 q = _recur4[0];
 y = _recur4[1];
 k = _recur4[2];

 y = (0, _jvs.jvs)(y, x) * q;
 return sign * y;
 }

 if (k <= 30.0) {
 k = 2.0;
 } else if (k < 90.0) {
 k = 3 * k / 4;
 }
 if (an > k + 3.0) {
 if (n < 0.0) k = -k;
 q = n - Math.floor(n);
 k = Math.floor(k) + q;
 if (n > 0.0) {
 var _recur5 = (0, _recur9.recur)(n, x, k, 1);

 var _recur6 = _slicedToArray(_recur5, 3);

 q = _recur6[0];
 n = _recur6[1];
 k = _recur6[2];
 } else {
 t = k;
 k = n;

 var _recur7 = (0, _recur9.recur)(t, x, k, 1);

 var _recur8 = _slicedToArray(_recur7, 3);

 q = _recur8[0];
 t = _recur8[1];
 k = _recur8[2];

 k = t;
 }
 if (q === 0.0) {
 y = 0.0;
 return sign * y;
 }
 } else {
 k = n;
 q = 1.0;
 }

 /* boundary between convergence of
 * power series and Hankel expansion
 */
 y = Math.abs(k);
 if (y < 26.0) {
 t = (0.0083 * y + 0.09) * y + 12.9;
 } else {
 t = 0.9 * y;
 }

 if (x > t) {
 y = (0, _hankel.hankel)(k, x);
 } else {
 y = (0, _jvs.jvs)(k, x);
 }
 if (n > 0.0) {
 y /= q;
 } else {
 y *= q;
 }
 } else {
 /* For large n, use the uniform expansion or the transitional expansion.
 * But if x is of the order of n**2, these may blow up, whereas the
 * Hankel expansion will then work.
 */
 if (n < 0.0) {
 // mtherr("Jv", TLOSS);
 y = NaN;
 return sign * y;
 }
 t = x / n;
 t /= n;
 if (t > 0.3) {
 y = (0, _hankel.hankel)(n, x);
 } else {
 y = (0, _jnx.jnx)(n, x);
 }
 }
 return sign * y;
}

exports.jv = jv;
},{"./constants.js":51,"./gamma.js":55,"./j0.js":67,"./j1.js":68,"./jv/hankel.js":70,"./jv/jnx.js":72,"./jv/jvs.js":73,"./jv/recur.js":74}],70:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.hankel = undefined;

var _constants = require('../constants.js');

var constants = _interopRequireWildcard(_constants);

function _interopRequireWildcard(obj) { if (obj && obj.__esModule) { return obj; } else { var newObj = {}; if (obj != null) { for (var key in obj) { if (Object.prototype.hasOwnProperty.call(obj, key)) newObj[key] = obj[key]; } } newObj.default = obj; return newObj; } }

/* Hankel's asymptotic expansion
* for large x.
* AMS55 #9.2.5.
*/

function hankel(n, x) {
 var t = void 0,
 u = void 0,
 z = void 0,
 k = void 0,
 sign = void 0,
 conv = void 0,
 p = void 0,
 q = void 0,
 j = void 0,
 m = void 0,
 pp = void 0,
 qq = void 0,
 flag = void 0;
 m = 4.0 * n * n;
 j = 1.0;
 z = 8.0 * x;
 k = 1.0;
 p = 1.0;
 u = (m - 1.0) / z;
 q = u;
 sign = 1.0;
 conv = 1.0;
 flag = 0;
 t = 1.0;
 pp = 1.0e38;
 qq = 1.0e38;

 while (t > constants.MACHEP) {
 k += 2.0;
 j += 1.0;
 sign = -sign;
 u *= (m - k * k) / (j * z);
 p += sign * u;
 k += 2.0;
 j += 1.0;
 u *= (m - k * k) / (j * z);
 q += sign * u;
 t = Math.abs(u / p);
 if (t < conv) {
 conv = t;
 qq = q;
 pp = p;
 flag = 1;
 }
 /* stop if the terms start getting larger */
 if (flag !== 0 && t > conv) break;
 }

 u = x - (0.5 * n + 0.25) * Math.PI;
 t = Math.sqrt(2.0 / (Math.PI * x)) * (pp * Math.cos(u) - qq * Math.sin(u));
 return t;
}

exports.hankel = hankel;
},{"../constants.js":51}],71:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.jnt = undefined;

var _slicedToArray = function () { function sliceIterator(arr, i) { var _arr = []; var _n = true; var _d = false; var _e = undefined; try { for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) { _arr.push(_s.value); if (i && _arr.length === i) break; } } catch (err) { _d = true; _e = err; } finally { try { if (!_n && _i["return"]) _i["return"](); } finally { if (_d) throw _e; } } return _arr; } return function (arr, i) { if (Array.isArray(arr)) { return arr; } else if (Symbol.iterator in Object(arr)) { return sliceIterator(arr, i); } else { throw new TypeError("Invalid attempt to destructure non-iterable instance"); } }; }();

var _airy3 = require('../airy.js');

var _polevl = require('../polevl.js');

/* Asymptotic expansion for transition region,
* n large and x close to n.
* AMS55 #9.3.23.
*/
var PF2 = [-9.0000000000000000000e-2, 8.5714285714285714286e-2];

var PF3 = [1.3671428571428571429e-1, -5.4920634920634920635e-2, -4.4444444444444444444e-3];

var PF4 = [1.3500000000000000000e-3, -1.6036054421768707483e-1, 4.2590187590187590188e-2, 2.7330447330447330447e-3];

var PG1 = [-2.4285714285714285714e-1, 1.4285714285714285714e-2];

var PG2 = [-9.0000000000000000000e-3, 1.9396825396825396825e-1, -1.1746031746031746032e-2];

var PG3 = [1.9607142857142857143e-2, -1.5983694083694083694e-1, 6.3838383838383838384e-3];

function jnt(n, x) {
 var z = void 0,
 zz = void 0,
 z3 = void 0,
 cbn = void 0,
 n23 = void 0,
 cbtwo = void 0,
 ai = void 0,
 aip = void 0,
 bi = void 0,
 bip = void 0,
 nk = void 0,
 fk = void 0,
 gk = void 0,
 pp = void 0,
 qq = void 0,
 k = void 0;

 var fArr = new Float64Array(5);
 var gArr = new Float64Array(4);

 cbn = Math.cbrt(n);
 z = (x - n) / cbn;
 cbtwo = Math.cbrt(2.0);

 /* Airy function */
 zz = -cbtwo * z;

 /* polynomials in expansion */
 var _airy = (0, _airy3.airy)(zz, ai, aip, bi, bip);

 var _airy2 = _slicedToArray(_airy, 4);

 ai = _airy2[0];
 aip = _airy2[1];
 bi = _airy2[2];
 bip = _airy2[3];
 zz = z * z;
 z3 = zz * z;
 fArr[0] = 1.0;
 fArr[1] = -z / 5.0;
 fArr[2] = (0, _polevl.polevl)(z3, PF2, 1) * zz;
 fArr[3] = (0, _polevl.polevl)(z3, PF3, 2);
 fArr[4] = (0, _polevl.polevl)(z3, PF4, 3) * z;
 gArr[0] = 0.3 * zz;
 gArr[1] = (0, _polevl.polevl)(z3, PG1, 1);
 gArr[2] = (0, _polevl.polevl)(z3, PG2, 2) * z;
 gArr[3] = (0, _polevl.polevl)(z3, PG3, 2) * zz;

 pp = 0.0;
 qq = 0.0;
 nk = 1.0;
 n23 = Math.cbrt(n * n);

 for (k = 0; k <= 4; k++) {
 fk = fArr[k] * nk;
 pp += fk;
 if (k !== 4) {
 gk = gArr[k] * nk;
 qq += gk;
 }

 nk /= n23;
 }

 fk = cbtwo * ai * pp / cbn + Math.cbrt(4.0) * aip * qq / n;
 return fk;
}

exports.jnt = jnt;
},{"../airy.js":49,"../polevl.js":78}],72:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.jnx = undefined;

var _slicedToArray = function () { function sliceIterator(arr, i) { var _arr = []; var _n = true; var _d = false; var _e = undefined; try { for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) { _arr.push(_s.value); if (i && _arr.length === i) break; } } catch (err) { _d = true; _e = err; } finally { try { if (!_n && _i["return"]) _i["return"](); } finally { if (_d) throw _e; } } return _arr; } return function (arr, i) { if (Array.isArray(arr)) { return arr; } else if (Symbol.iterator in Object(arr)) { return sliceIterator(arr, i); } else { throw new TypeError("Invalid attempt to destructure non-iterable instance"); } }; }();

var _airy3 = require('../airy.js');

var _polevl = require('../polevl.js');

var _jnt = require('./jnt.js');

var _constants = require('../constants.js');

var constants = _interopRequireWildcard(_constants);

function _interopRequireWildcard(obj) { if (obj && obj.__esModule) { return obj; } else { var newObj = {}; if (obj != null) { for (var key in obj) { if (Object.prototype.hasOwnProperty.call(obj, key)) newObj[key] = obj[key]; } } newObj.default = obj; return newObj; } }

/* Asymptotic expansion for large n.
* AMS55 #9.3.35.
*/
var lambda = [1.0, 1.041666666666666666666667E-1, 8.355034722222222222222222E-2, 1.282265745563271604938272E-1, 2.918490264641404642489712E-1, 8.816272674437576524187671E-1, 3.321408281862767544702647E+0, 1.499576298686255465867237E+1, 7.892301301158651813848139E+1, 4.744515388682643231611949E+2, 3.207490090890661934704328E+3];

var mu = [1.0, -1.458333333333333333333333E-1, -9.874131944444444444444444E-2, -1.433120539158950617283951E-1, -3.172272026784135480967078E-1, -9.424291479571202491373028E-1, -3.511203040826354261542798E+0, -1.572726362036804512982712E+1, -8.228143909718594444224656E+1, -4.923553705236705240352022E+2, -3.316218568547972508762102E+3];

var P1 = [-2.083333333333333333333333E-1, 1.250000000000000000000000E-1];

var P2 = [3.342013888888888888888889E-1, -4.010416666666666666666667E-1, 7.031250000000000000000000E-2];

var P3 = [-1.025812596450617283950617E+0, 1.846462673611111111111111E+0, -8.912109375000000000000000E-1, 7.324218750000000000000000E-2];

var P4 = [4.669584423426247427983539E+0, -1.120700261622299382716049E+1, 8.789123535156250000000000E+0, -2.364086914062500000000000E+0, 1.121520996093750000000000E-1];

var P5 = [-2.8212072558200244877E1, 8.4636217674600734632E1, -9.1818241543240017361E1, 4.2534998745388454861E1, -7.3687943594796316964E0, 2.27108001708984375E-1];

var P6 = [2.1257013003921712286E2, -7.6525246814118164230E2, 1.0599904525279998779E3, -6.9957962737613254123E2, 2.1819051174421159048E2, -2.6491430486951555525E1, 5.7250142097473144531E-1];

var P7 = [-1.9194576623184069963E3, 8.0617221817373093845E3, -1.3586550006434137439E4, 1.1655393336864533248E4, -5.3056469786134031084E3, 1.2009029132163524628E3, -1.0809091978839465550E2, 1.7277275025844573975E0];

function jnx(n, x) {
 var zeta = void 0,
 sqz = void 0,
 zz = void 0,
 zp = void 0,
 np = void 0,
 cbn = void 0,
 n23 = void 0,
 t = void 0,
 z = void 0,
 sz = void 0,
 pp = void 0,
 qq = void 0,
 z32i = void 0,
 zzi = void 0,
 ak = void 0,
 bk = void 0,
 akl = void 0,
 bkl = void 0,
 sign = void 0,
 doa = void 0,
 dob = void 0,
 nflg = void 0,
 k = void 0,
 s = void 0,
 tk = void 0,
 tkp1 = void 0,
 m = void 0,
 ai = void 0,
 aip = void 0,
 bi = void 0,
 bip = void 0,
 uArr = void 0;
 uArr = new Float64Array(8);

 /* Test for x very close to n. Use expansion for transition region if so. */
 cbn = Math.cbrt(n);
 z = (x - n) / cbn;
 if (Math.abs(z) <= 0.7) return (0, _jnt.jnt)(n, x);

 z = x / n;
 zz = 1.0 - z * z;
 if (zz === 0.0) return 0.0;

 if (zz > 0.0) {
 sz = Math.sqrt(zz);
 t = 1.5 * (Math.log((1.0 + sz) / z) - sz); /* zeta ** 3/2 */
 zeta = Math.cbrt(t * t);
 nflg = 1;
 } else {
 sz = Math.sqrt(-zz);
 t = 1.5 * (sz - Math.acos(1.0 / z));
 zeta = -Math.cbrt(t * t);
 nflg = -1;
 }
 z32i = Math.abs(1.0 / t);
 sqz = Math.cbrt(t);

 /* Airy function */
 n23 = Math.cbrt(n * n);
 t = n23 * zeta;

 /* polynomials in expansion */
 var _airy = (0, _airy3.airy)(t, ai, aip, bi, bip);

 var _airy2 = _slicedToArray(_airy, 4);

 ai = _airy2[0];
 aip = _airy2[1];
 bi = _airy2[2];
 bip = _airy2[3];
 uArr[0] = 1.0;
 zzi = 1.0 / zz;
 uArr[1] = (0, _polevl.polevl)(zzi, P1, 1) / sz;
 uArr[2] = (0, _polevl.polevl)(zzi, P2, 2) / zz;
 uArr[3] = (0, _polevl.polevl)(zzi, P3, 3) / (sz * zz);
 pp = zz * zz;
 uArr[4] = (0, _polevl.polevl)(zzi, P4, 4) / pp;
 uArr[5] = (0, _polevl.polevl)(zzi, P5, 5) / (pp * sz);
 pp *= zz;
 uArr[6] = (0, _polevl.polevl)(zzi, P6, 6) / pp;
 uArr[7] = (0, _polevl.polevl)(zzi, P7, 7) / (pp * sz);

 pp = 0.0;
 qq = 0.0;
 np = 1.0;
 /* flags to stop when terms get larger */
 doa = 1;
 dob = 1;
 akl = Infinity;
 bkl = Infinity;

 for (k = 0; k <= 3; k++) {
 tk = 2 * k;
 tkp1 = tk + 1;
 zp = 1.0;
 ak = 0.0;
 bk = 0.0;
 for (s = 0; s <= tk; s++) {
 if (doa) {
 if ((s & 3) > 1) sign = nflg;else sign = 1;
 ak += sign * mu[s] * zp * uArr[tk - s];
 }

 if (dob) {
 m = tkp1 - s;
 if ((m + 1 & 3) > 1) sign = nflg;else sign = 1;
 bk += sign * lambda[s] * zp * uArr[m];
 }
 zp *= z32i;
 }

 if (doa) {
 ak *= np;
 t = Math.abs(ak);
 if (t < akl) {
 akl = t;
 pp += ak;
 } else {
 doa = 0;
 }
 }

 if (dob) {
 bk += lambda[tkp1] * zp * uArr[0];
 bk *= -np / sqz;
 t = Math.abs(bk);
 if (t < bkl) {
 bkl = t;
 qq += bk;
 } else {
 dob = 0;
 }
 }

 if (np < constants.MACHEP) break;
 np /= n * n;
 }

 /* normalizing factor ( 4*zeta/(1 - z**2) )**1/4 */
 t = 4.0 * zeta / zz;
 t = Math.sqrt(Math.sqrt(t));

 t *= ai * pp / Math.cbrt(n) + aip * qq / (n23 * n);
 return t;
}

exports.jnx = jnx;
},{"../airy.js":49,"../constants.js":51,"../polevl.js":78,"./jnt.js":71}],73:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.jvs = undefined;

var _slicedToArray = function () { function sliceIterator(arr, i) { var _arr = []; var _n = true; var _d = false; var _e = undefined; try { for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) { _arr.push(_s.value); if (i && _arr.length === i) break; } } catch (err) { _d = true; _e = err; } finally { try { if (!_n && _i["return"]) _i["return"](); } finally { if (_d) throw _e; } } return _arr; } return function (arr, i) { if (Array.isArray(arr)) { return arr; } else if (Symbol.iterator in Object(arr)) { return sliceIterator(arr, i); } else { throw new TypeError("Invalid attempt to destructure non-iterable instance"); } }; }();

var _constants = require('../constants.js');

var constants = _interopRequireWildcard(_constants);

var _pythonHelpers = require('../../../utils/pythonHelpers.js');

var py = _interopRequireWildcard(_pythonHelpers);

var _gamma = require('../gamma.js');

var _lgam = require('../gamma/lgam.js');

function _interopRequireWildcard(obj) { if (obj && obj.__esModule) { return obj; } else { var newObj = {}; if (obj != null) { for (var key in obj) { if (Object.prototype.hasOwnProperty.call(obj, key)) newObj[key] = obj[key]; } } newObj.default = obj; return newObj; } }

/* Ascending power series for Jv(x).
* AMS55 #9.1.10.
*/
function jvs(n, x) {
 var t = void 0,
 u = void 0,
 y = void 0,
 z = void 0,
 k = void 0,
 ex = void 0,
 sgngam = void 0,
 res = void 0;

 z = -x * x / 4.0;
 u = 1.0;
 y = u;
 k = 1.0;
 t = 1.0;

 while (t > constants.MACHEP) {
 u *= z / (k * (n + k));
 y += u;
 k += 1.0;
 if (y !== 0) t = Math.abs(u / y);
 }

 var _py$frexp = py.frexp(0.5 * x);

 var _py$frexp2 = _slicedToArray(_py$frexp, 2);

 t = _py$frexp2[0];
 ex = _py$frexp2[1];

 ex = ex * n;
 if (ex > -1023 && ex < 1023 && n > 0.0 && n < constants.MAXGAM - 1.0) {
 t = Math.pow(0.5 * x, n) / (0, _gamma.gamma)(n + 1.0);
 y *= t;
 } else {
 var _lgamSgn = (0, _lgam.lgamSgn)(n + 1.0);

 var _lgamSgn2 = _slicedToArray(_lgamSgn, 2);

 res = _lgamSgn2[0];
 sgngam = _lgamSgn2[1];

 t = n * Math.log(0.5 * x) - res;
 if (y < 0) {
 sgngam = -sgngam;
 y = -y;
 }
 t += Math.log(y);
 if (t < -constants.MAXLOG) {
 return 0.0;
 }
 if (t > constants.MAXLOG) {
 // mtherr("Jv", OVERFLOW);
 return Infinity;
 }
 y = sgngam * Math.exp(t);
 }
 return y;
}

exports.jvs = jvs;
},{"../../../utils/pythonHelpers.js":94,"../constants.js":51,"../gamma.js":55,"../gamma/lgam.js":56}],74:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.recur = undefined;

var _constants = require('../constants.js');

var constants = _interopRequireWildcard(_constants);

function _interopRequireWildcard(obj) { if (obj && obj.__esModule) { return obj; } else { var newObj = {}; if (obj != null) { for (var key in obj) { if (Object.prototype.hasOwnProperty.call(obj, key)) newObj[key] = obj[key]; } } newObj.default = obj; return newObj; } }

/* Reduce the order by backward recurrence.
* AMS55 #9.1.27 and 9.1.73.
*/
function recur(n, x, newn, cancel) {
 var pkm2 = void 0,
 pkm1 = void 0,
 pk = void 0,
 qkm2 = void 0,
 qkm1 = void 0,
 k = void 0,
 ans = void 0,
 qk = void 0,
 xk = void 0,
 yk = void 0,
 r = void 0,
 t = void 0,
 kf = void 0,
 nflag = void 0,
 ctr = void 0,
 miniter = void 0,
 maxiter = void 0;

 /* Continued fraction for Jn(x)/Jn-1(x)
 * AMS 9.1.73
 *
 * x -x^2 -x^2
 * ------ --------- --------- ...
 * 2 n + 2(n+1) + 2(n+2) +
 *
 * Compute it with the simplest possible algorithm.
 *
 * This continued fraction starts to converge when (|n| + m) > |x|.
 * Hence, at least |x|-|n| iterations are necessary before convergence is
 * achieved. There is a hard limit set below, m <= 30000, which is chosen
 * so that no branch in `jv` requires more iterations to converge.
 * The exact maximum number is (500/3.6)^2 - 500 ~ 19000
 */
 var BIG = 1.44115188075855872E+17;
 maxiter = 22000;
 miniter = Math.abs(x) - Math.abs(n);
 if (miniter < 1) miniter = 1;
 if (n < 0.0) nflag = 1;else nflag = 0;

 var goToLabel = 'fstart';
 mainExecutionLoop: while (true) {
 innerSwitch: switch (goToLabel) {
 case 'fstart':
 pkm2 = 0.0;
 qkm2 = 1.0;
 pkm1 = x;
 qkm1 = 2 * n;
 xk = -x * x;
 yk = qkm1;
 ans = 0.0; /* ans=0.0 ensures that t=1.0 in the first iteration */
 ctr = 0;
 do {
 yk += 2.0;
 pk = pkm1 * yk + pkm2 * xk;
 qk = qkm1 * yk + qkm2 * xk;
 pkm2 = pkm1;
 pkm1 = pk;
 qkm2 = qkm1;
 qkm1 = qk;

 /* check convergence */
 if (qk !== 0 && ctr > miniter) r = pk / qk;else r = 0.0;

 if (r !== 0) {
 t = Math.abs((ans - r) / r);
 ans = r;
 } else {
 t = 1.0;
 }

 if (++ctr > maxiter) {
 // mtherr("jv", UNDERFLOW);
 goToLabel = 'done';break innerSwitch;
 }
 if (t < constants.MACHEP) {
 goToLabel = 'done';break innerSwitch;
 }

 /* renormalize coefficients */
 if (Math.abs(pk) > BIG) {
 pkm2 /= BIG;
 pkm1 /= BIG;
 qkm2 /= BIG;
 qkm1 /= BIG;
 }
 } while (t > constants.MACHEP);

 case 'done':
 if (ans === 0) ans = 1.0;

 /* Change n to n-1 if n < 0 and the continued fraction is small */
 if (nflag > 0) {
 if (Math.abs(ans) < 0.125) {
 nflag = -1;
 n--;
 goToLabel = 'fstart';break;
 }
 }

 kf = newn;

 /* backward recurrence
 * 2k
 * J (x) = --- J (x) - J (x)
 * k-1 x k k+1
 */

 pk = 1.0;
 pkm1 = 1.0 / ans;
 k = n - 1.0;
 r = 2 * k;
 do {
 pkm2 = (pkm1 * r - pk * x) / x;
 /* pkp1 = pk; */
 pk = pkm1;
 pkm1 = pkm2;
 r -= 2.0;
 /*
 * t = Math.abs(pkp1) + Math.abs(pk);
 * if( (k > (kf + 2.5)) && (Math.abs(pkm1) < 0.25*t) )
 * {
 * k -= 1.0;
 * t = x*x;
 * pkm2 = ( (r*(r+2.0)-t)*pk - r*x*pkp1 )/t;
 * pkp1 = pk;
 * pk = pkm1;
 * pkm1 = pkm2;
 * r -= 2.0;
 * }
 */
 k -= 1.0;
 } while (k > kf + 0.5);

 /* Take the larger of the last two iterates
 * on the theory that it may have less cancellation error.
 */

 if (cancel) {
 if (kf >= 0.0 && Math.abs(pk) > Math.abs(pkm1)) {
 k += 1.0;
 pkm2 = pk;
 }
 }
 newn = k;
 default:
 break mainExecutionLoop;
 }
 }
 return [pkm2, n, newn];
} /* eslint-disable no-labels */
/* eslint-disable no-fallthrough */
exports.recur = recur;
},{"../constants.js":51}],75:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.k0e = exports.k0 = undefined;

var _chbevl = require('./chbevl.js');

var _i = require('./i0.js');

/**
 * @file k0.js Modified Bessel function, third kind, order zero
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, k0();
 *
 * y = k0( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns modified Bessel function of the third kind
 * of order zero of the argument.
 *
 * The range is partitioned into the two intervals [0,8] and
 * (8, infinity). Chebyshev polynomial expansions are employed
 * in each interval.
 *
 *
 *
 * ACCURACY:
 *
 * Tested at 2000 random points between 0 and 8. Peak absolute
 * error (relative when K0 > 1) was 1.46e-14; rms, 4.26e-15.
 * Relative error:
 * arithmetic domain # trials peak rms
 * IEEE 0, 30 30000 1.2e-15 1.6e-16
 *
 * ERROR MESSAGES:
 *
 * message condition value returned
 * K0 domain x <= 0 NPY_INFINITY
 *
 *
 * Modified Bessel function, third kind, order zero,
 * exponentially scaled
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, k0e();
 *
 * y = k0e( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns exponentially scaled modified Bessel function
 * of the third kind of order zero of the argument.
 *
 *
 *
 * ACCURACY:
 *
 * Relative error:
 * arithmetic domain # trials peak rms
 * IEEE 0, 30 30000 1.4e-15 1.4e-16
 * See k0().
 *
 * Cephes Math Library Release 2.8: June, 2000
 * Copyright 1984, 1987, 2000 by Stephen L. Moshier
 * Ported to ECMAScript 2018
 * Copyright (c) 2018, Kings Distributed Systems
 *
 * @author KC Erb, kc@kcerb.com
 * @date April 2018
 */
var A = [1.37446543561352307156E-16, 4.25981614279661018399E-14, 1.03496952576338420167E-11, 1.90451637722020886025E-9, 2.53479107902614945675E-7, 2.28621210311945178607E-5, 1.26461541144692592338E-3, 3.59799365153615016266E-2, 3.44289899924628486886E-1, -5.35327393233902768720E-1];

/* Chebyshev coefficients for exp(x) sqrt(x) K0(x)
 * in the inverted interval [2,infinity].
 *
 * lim(x->inf){ exp(x) sqrt(x) K0(x) } = sqrt(pi/2).
 */
var B = [5.30043377268626276149E-18, -1.64758043015242134646E-17, 5.21039150503902756861E-17, -1.67823109680541210385E-16, 5.51205597852431940784E-16, -1.84859337734377901440E-15, 6.34007647740507060557E-15, -2.22751332699166985548E-14, 8.03289077536357521100E-14, -2.98009692317273043925E-13, 1.14034058820847496303E-12, -4.51459788337394416547E-12, 1.85594911495471785253E-11, -7.95748924447710747776E-11, 3.57739728140030116597E-10, -1.69753450938905987466E-9, 8.57403401741422608519E-9, -4.66048989768794782956E-8, 2.76681363944501510342E-7, -1.83175552271911948767E-6, 1.39498137188764993662E-5, -1.28495495816278026384E-4, 1.56988388573005337491E-3, -3.14481013119645005427E-2, 2.44030308206595545468E0];

function k0(x) {
 var y = void 0,
 z = void 0;

 if (x === 0.0) {
 // mtherr("k0", SING);
 return Infinity;
 } else if (x < 0.0) {
 // mtherr("k0", DOMAIN);
 return NaN;
 }

 if (x <= 2.0) {
 y = x * x - 2.0;
 y = (0, _chbevl.chbevl)(y, A, 10) - Math.log(0.5 * x) * (0, _i.i0)(x);
 return y;
 }
 z = 8.0 / x - 2.0;
 y = Math.exp(-x) * (0, _chbevl.chbevl)(z, B, 25) / Math.sqrt(x);
 return y;
}

function k0e(x) {
 var y = void 0;

 if (x === 0.0) {
 // mtherr("k0e", SING);
 return Infinity;
 } else if (x < 0.0) {
 // mtherr("k0e", DOMAIN);
 return NaN;
 }

 if (x <= 2.0) {
 y = x * x - 2.0;
 y = (0, _chbevl.chbevl)(y, A, 10) - Math.log(0.5 * x) * (0, _i.i0)(x);
 return y * Math.exp(x);
 }

 y = (0, _chbevl.chbevl)(8.0 / x - 2.0, B, 25) / Math.sqrt(x);
 return y;
}

exports.k0 = k0;
exports.k0e = k0e;
},{"./chbevl.js":50,"./i0.js":62}],76:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.k1e = exports.k1 = undefined;

var _chbevl = require('./chbevl.js');

var _i = require('./i1.js');

/**
 * @file k1.js Modified Bessel function, third kind, order one
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, k1();
 *
 * y = k1( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Computes the modified Bessel function of the third kind
 * of order one of the argument.
 *
 * The range is partitioned into the two intervals [0,2] and
 * (2, infinity). Chebyshev polynomial expansions are employed
 * in each interval.
 *
 *
 *
 * ACCURACY:
 *
 * Relative error:
 * arithmetic domain # trials peak rms
 * IEEE 0, 30 30000 1.2e-15 1.6e-16
 *
 * ERROR MESSAGES:
 *
 * message condition value returned
 * k1 domain x <= 0 NPY_INFINITY
 *
 *
 * Modified Bessel function, third kind, order one,
 * exponentially scaled
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, k1e();
 *
 * y = k1e( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns exponentially scaled modified Bessel function
 * of the third kind of order one of the argument:
 *
 * k1e(x) = exp(x) * k1(x).
 *
 *
 *
 * ACCURACY:
 *
 * Relative error:
 * arithmetic domain # trials peak rms
 * IEEE 0, 30 30000 7.8e-16 1.2e-16
 * See k1().
 *
 * Cephes Math Library Release 2.8: June, 2000
 * Copyright 1984, 1987, 2000 by Stephen L. Moshier
 * Ported to ECMAScript 2018
 * Copyright (c) 2018, Kings Distributed Systems
 *
 * @author KC Erb, kc@kcerb.com
 * @date April 2018
 */
var A = [-7.02386347938628759343E-18, -2.42744985051936593393E-15, -6.66690169419932900609E-13, -1.41148839263352776110E-10, -2.21338763073472585583E-8, -2.43340614156596823496E-6, -1.73028895751305206302E-4, -6.97572385963986435018E-3, -1.22611180822657148235E-1, -3.53155960776544875667E-1, 1.52530022733894777053E0];

/* Chebyshev coefficients for exp(x) sqrt(x) K1(x)
* in the interval [2,infinity].
*
* lim(x->inf){ exp(x) sqrt(x) K1(x) } = sqrt(pi/2).
*/
var B = [-5.75674448366501715755E-18, 1.79405087314755922667E-17, -5.68946255844285935196E-17, 1.83809354436663880070E-16, -6.05704724837331885336E-16, 2.03870316562433424052E-15, -7.01983709041831346144E-15, 2.47715442448130437068E-14, -8.97670518232499435011E-14, 3.34841966607842919884E-13, -1.28917396095102890680E-12, 5.13963967348173025100E-12, -2.12996783842756842877E-11, 9.21831518760500529508E-11, -4.19035475934189648750E-10, 2.01504975519703286596E-9, -1.03457624656780970260E-8, 5.74108412545004946722E-8, -3.50196060308781257119E-7, 2.40648494783721712015E-6, -1.93619797416608296024E-5, 1.95215518471351631108E-4, -2.85781685962277938680E-3, 1.03923736576817238437E-1, 2.72062619048444266945E0];

function k1(x) {
 var y = void 0,
 z = void 0;

 if (x === 0.0) {
 // mtherr("k1", SING);
 return Infinity;
 } else if (x < 0.0) {
 // mtherr("k1", DOMAIN);
 return NaN;
 }
 z = 0.5 * x;

 if (x <= 2.0) {
 y = x * x - 2.0;
 y = Math.log(z) * (0, _i.i1)(x) + (0, _chbevl.chbevl)(y, A, 11) / x;
 return y;
 }

 return Math.exp(-x) * (0, _chbevl.chbevl)(8.0 / x - 2.0, B, 25) / Math.sqrt(x);
}

function k1e(x) {
 var y = void 0;

 if (x === 0.0) {
 // mtherr("k1e", SING);
 return Infinity;
 } else if (x < 0.0) {
 // mtherr("k1e", DOMAIN);
 return NaN;
 }

 if (x <= 2.0) {
 y = x * x - 2.0;
 y = Math.log(0.5 * x) * (0, _i.i1)(x) + (0, _chbevl.chbevl)(y, A, 11) / x;
 return y * Math.exp(x);
 }

 return (0, _chbevl.chbevl)(8.0 / x - 2.0, B, 25) / Math.sqrt(x);
}

exports.k1 = k1;
exports.k1e = k1e;
},{"./chbevl.js":50,"./i1.js":63}],77:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.kn = undefined;

var _constants = require('./constants.js');

var constants = _interopRequireWildcard(_constants);

function _interopRequireWildcard(obj) { if (obj && obj.__esModule) { return obj; } else { var newObj = {}; if (obj != null) { for (var key in obj) { if (Object.prototype.hasOwnProperty.call(obj, key)) newObj[key] = obj[key]; } } newObj.default = obj; return newObj; } }

var EUL = 5.772156649015328606065e-1; /**
 * @file kn.js Modified Bessel function, third kind, integer order
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, kn();
 * int n;
 *
 * y = kn( n, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns modified Bessel function of the third kind
 * of order n of the argument.
 *
 * The range is partitioned into the two intervals [0,9.55] and
 * (9.55, infinity). An ascending power series is used in the
 * low range, and an asymptotic expansion in the high range.
 *
 *
 *
 * ACCURACY:
 *
 * Relative error:
 * arithmetic domain # trials peak rms
 * IEEE 0,30 90000 1.8e-8 3.0e-10
 *
 * Error is high only near the crossover point x = 9.55
 * between the two expansions used.
 *
 * Cephes Math Library Release 2.8: June, 2000
 * Copyright 1984, 1987, 1988, 2000 by Stephen L. Moshier
 * Ported to ECMAScript 2018
 * Copyright (c) 2018, Kings Distributed Systems
 *
 * @author KC Erb, kc@kcerb.com
 * @date April 2018
 */

/*
 * Algorithm for Kn.
 * n-1
 * -n - (n-k-1)! 2 k
 * K (x) = 0.5 (x/2) > -------- (-x /4)
 * n - k!
 * k=0
 *
 * inf. 2 k
 * n n - (x /4)
 * + (-1) 0.5(x/2) > {p(k+1) + p(n+k+1) - 2log(x/2)} ---------
 * - k! (n+k)!
 * k=0
 *
 * where p(m) is the psi function: p(1) = -EUL and
 *
 * m-1
 * -
 * p(m) = -EUL + > 1/k
 * -
 * k=1
 *
 * For large x,
 * 2 2 2
 * u-1 (u-1 )(u-3 )
 * K (z) = sqrt(pi/2z) exp(-z) { 1 + ------- + ------------ + ...}
 * v 1 2
 * 1! (8z) 2! (8z)
 * asymptotically, where
 *
 * 2
 * u = 4 v .
 *
 */

var MAXFAC = 31;

function kn(nn, x) {
 var k = void 0,
 kf = void 0,
 nk1f = void 0,
 zn = void 0,
 t = void 0,
 s = void 0,
 z0 = void 0,
 z = void 0,
 ans = void 0,
 fn = void 0,
 pn = void 0,
 pk = void 0,
 zmn = void 0,
 tlg = void 0,
 tox = void 0,
 i = void 0,
 n = void 0;

 if (nn < 0) n = -nn;else n = nn;

 if (n > MAXFAC) {
 // mtherr("kn", OVERFLOW);
 return Infinity;
 }

 if (x <= 0.0) {
 if (x < 0.0) {
 // mtherr("kn", DOMAIN);
 return NaN;
 } else {
 // mtherr("kn", SING);
 return Infinity;
 }
 }

 if (x > 9.55) return asymp(x, n);

 ans = 0.0;
 z0 = 0.25 * x * x;
 fn = 1.0;
 pn = 0.0;
 zmn = 1.0;
 tox = 2.0 / x;

 if (n > 0) {
 /* compute factorial of n and psi(n) */
 pn = -EUL;
 k = 1.0;
 for (i = 1; i < n; i++) {
 pn += 1.0 / k;
 k += 1.0;
 fn *= k;
 }

 zmn = tox;

 if (n === 1) {
 ans = 1.0 / x;
 } else {
 nk1f = fn / n;
 kf = 1.0;
 s = nk1f;
 z = -z0;
 zn = 1.0;
 for (i = 1; i < n; i++) {
 nk1f = nk1f / (n - i);
 kf = kf * i;
 zn *= z;
 t = nk1f * zn / kf;
 s += t;
 if (constants.DBL_MAX - Math.abs(t) < Math.abs(s)) {
 // mtherr("kn", OVERFLOW);
 return Infinity;
 }
 if (tox > 1.0 && constants.DBL_MAX / tox < zmn) {
 // mtherr("kn", OVERFLOW);
 return Infinity;
 }
 zmn *= tox;
 }
 s *= 0.5;
 t = Math.abs(s);
 if (zmn > 1.0 && constants.DBL_MAX / zmn < t) {
 // mtherr("kn", OVERFLOW);
 return Infinity;
 }
 if (t > 1.0 && constants.DBL_MAX / t < zmn) {
 // mtherr("kn", OVERFLOW);
 return Infinity;
 }
 ans = s * zmn;
 }
 }

 tlg = 2.0 * Math.log(0.5 * x);
 pk = -EUL;
 if (n === 0) {
 pn = pk;
 t = 1.0;
 } else {
 pn = pn + 1.0 / n;
 t = 1.0 / fn;
 }
 s = (pk + pn - tlg) * t;
 k = 1.0;
 do {
 t *= z0 / (k * (k + n));
 pk += 1.0 / k;
 pn += 1.0 / (k + n);
 s += (pk + pn - tlg) * t;
 k += 1.0;
 } while (Math.abs(t / s) > constants.MACHEP);

 s = 0.5 * s / zmn;
 if (n & 1) s = -s;
 ans += s;

 return ans;
}

/* Asymptotic expansion for Kn(x) */
/* Converges to 1.4e-17 for x > 18.4 */
function asymp(x, n) {
 var nkf = void 0,
 k = void 0,
 pn = void 0,
 pk = void 0,
 z0 = void 0,
 fn = void 0,
 t = void 0,
 s = void 0,
 i = void 0,
 z = void 0,
 nk1f = void 0,
 ans = void 0;
 if (x > constants.MAXLOG) {
 // mtherr("kn", UNDERFLOW);
 return 0.0;
 }
 k = n;
 pn = 4.0 * k * k;
 pk = 1.0;
 z0 = 8.0 * x;
 fn = 1.0;
 t = 1.0;
 s = t;
 nkf = Infinity;
 i = 0;
 do {
 z = pn - pk * pk;
 t = t * z / (fn * z0);
 nk1f = Math.abs(t);
 if (i >= n && nk1f > nkf) break;
 nkf = nk1f;
 s += t;
 fn += 1.0;
 pk += 2.0;
 i += 1;
 } while (Math.abs(t / s) > constants.MACHEP);

 ans = Math.exp(-x) * Math.sqrt(Math.PI / (2.0 * x)) * s;
 return ans;
}

exports.kn = kn;
},{"./constants.js":51}],78:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
/**
 * @file polevl.js Evaluate polynomial
 *
 *
 *
 * SYNOPSIS:
 *
 * int N;
 * double x, y, coef[N+1], polevl[];
 *
 * y = polevl( x, coef, N );
 *
 *
 *
 * DESCRIPTION:
 *
 * Evaluates polynomial of degree N:
 *
 * 2 N
 * y = C + C x + C x +...+ C x
 * 0 1 2 N
 *
 * Coefficients are stored in reverse order:
 *
 * coef[0] = C , ..., coef[N] = C .
 * N 0
 *
 * The function p1evl() assumes that coef[N] = 1.0 and is
 * omitted from the array. Its calling arguments are
 * otherwise the same as polevl().
 *
 *
 * SPEED:
 *
 * In the interest of speed, there are no checks for out
 * of bounds arithmetic. This routine is used by most of
 * the functions in the library. Depending on available
 * equipment features, the user may wish to rewrite the
 * program in microcode or assembly language.
 *
 * Cephes Math Library Release 2.1: December, 1988
 * Copyright 1984, 1987, 1988 by Stephen L. Moshier
 * Ported to ECMAScript 2018
 * Copyright (c) 2018, Kings Distributed Systems
 *
 * @author KC Erb, kc@kcerb.com
 * @date April 2018
 */
function polevl(x, coef, n) {
 var ans = void 0,
 i = void 0;
 ans = coef[0];
 for (i = 1; i <= n; i++) {
 ans = ans * x + coef[i];
 }

 return ans;
}

function p1evl(x, coef, n) {
 var ans = void 0,
 i = void 0;

 ans = x + coef[0];

 for (i = 1; i < n; i++) {
 ans = ans * x + coef[i];
 }

 return ans;
};

// TODO: ratevl

exports.polevl = polevl;
exports.p1evl = p1evl;
},{}],79:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.psi = undefined;

var _polevl = require('./polevl.js');

// Psi (digamma) function
function psi(x) {
 var q = void 0,
 p = void 0;
 var A = new Float64Array([8.33333333333333333333E-2, -8.33333333333333333333E-3, 3.96825396825396825397E-3, -4.16666666666666666667E-3, 7.57575757575757575758E-3, -2.10927960927960927961E-2, 8.33333333333333333333E-2]);

 var negative = 0;
 var nz = 0.0;
 var y = void 0;

 if (x <= 0.0) {
 negative = 1;
 q = x;
 p = Math.floor(q);
 if (p === q) {
 // mtherr("psi", SING);
 return Infinity;
 }
 // Remove the zeros of tan(PI x)
 // by subtracting the nearest integer from x
 nz = q - p;
 if (nz !== 0.5) {
 if (nz > 0.5) {
 p += 1.0;
 nz = q - p;
 }
 nz = Math.PI / Math.tan(Math.PI * nz);
 } else {
 nz = 0.0;
 }
 x = 1.0 - x;
 }

 // check for positive integer up to 10
 if (x <= 10.0 && x === Math.floor(x)) {
 y = 0.0;
 var n = x;
 for (var i = 1; i < n; i++) {
 var _w = i;
 y += 1.0 / _w;
 }
 y -= Math.E;
 return done(negative, y, nz);
 }

 var s = x;
 var w = 0.0;
 while (s < 10.0) {
 w += 1.0 / s;
 s += 1.0;
 }

 if (s < 1.0e17) {
 var z = 1.0 / (s * s);
 y = z * (0, _polevl.polevl)(z, A, 6);
 } else {
 y = 0.0;
 }

 y = Math.log(s) - 0.5 / s - y - w;

 return done(negative, y, nz);
} /**
 * @file psi.js Psi (digamma) function
 *
 *
 * SYNOPSIS:
 *
 * double x, y, psi();
 *
 * y = psi( x );
 *
 *
 * DESCRIPTION:
 *
 * d -
 * psi(x) = -- ln | (x)
 * dx
 *
 * is the logarithmic derivative of the gamma function.
 * For integer x,
 * n-1
 * -
 * psi(n) = -EUL + > 1/k.
 * -
 * k=1
 *
 * This formula is used for 0 < n <= 10. If x is negative, it
 * is transformed to a positive argument by the reflection
 * formula psi(1-x) = psi(x) + pi cot(pi x).
 * For general positive x, the argument is made greater than 10
 * using the recurrence psi(x+1) = psi(x) + 1/x.
 * Then the following asymptotic expansion is applied:
 *
 * inf. B
 * - 2k
 * psi(x) = log(x) - 1/2x - > -------
 * - 2k
 * k=1 2k x
 *
 * where the B2k are Bernoulli numbers.
 *
 * ACCURACY:
 * Relative error (except absolute when |psi| < 1):
 * arithmetic domain # trials peak rms
 * IEEE 0,30 30000 1.3e-15 1.4e-16
 * IEEE -30,0 40000 1.5e-15 2.2e-16
 *
 * ERROR MESSAGES:
 * message condition value returned
 * psi singularity x integer <=0 NPY_INFINITY
 *
 * Cephes Math Library Release 2.8: June, 2000
 * Copyright 1984, 1987, 1992, 2000 by Stephen L. Moshier
 *
 * Code for the rational approximation on [1, 2] is:
 *
 * (C) Copyright John Maddock 2006.
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 * Ported to ECMAScript 2018
 * Copyright (c) 2018, Kings Distributed Systems
 *
 * @author KC Erb, kc@kcerb.com
 * @date April 2018
 */
;

function done(negative, y, nz) {
 if (negative) {
 y -= nz;
 }
 return y;
};

exports.psi = psi;
},{"./polevl.js":78}],80:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.yn = undefined;

var _j = require('./j0.js');

var _j2 = require('./j1.js');

/**
 * @file yn.js Bessel function of second kind of integer order
 *
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, yn();
 * int n;
 *
 * y = yn( n, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of order n, where n is a
 * (possibly negative) integer.
 *
 * The function is evaluated by forward recurrence on
 * n, starting with values computed by the routines
 * y0() and y1().
 *
 * If n = 0 or 1 the routine for y0 or y1 is called
 * directly.
 *
 *
 *
 * ACCURACY:
 *
 *
 * Absolute error, except relative
 * when y > 1:
 * arithmetic domain # trials peak rms
 * IEEE 0, 30 30000 3.4e-15 4.3e-16
 *
 *
 * ERROR MESSAGES:
 *
 * message condition value returned
 * yn singularity x = 0 NPY_INFINITY
 * yn overflow NPY_INFINITY
 *
 * Spot checked against tables for x, n between 0 and 100.
 *
 *
 * Cephes Math Library Release 2.8: June, 2000
 * Copyright 1984, 1987, 2000 by Stephen L. Moshier
 * Ported to ECMAScript 2018
 * Copyright (c) 2018, Kings Distributed Systems
 *
 * @author KC Erb, kc@kcerb.com
 * @date April 2018
 */

function yn(n, x) {
 var an = void 0,
 anm1 = void 0,
 anm2 = void 0,
 r = void 0,
 k = void 0,
 sign = void 0;

 if (n < 0) {
 n = -n;
 /* -1**n */
 if ((n & 1) === 0) {
 sign = 1;
 } else {
 sign = -1;
 }
 } else {
 sign = 1;
 }

 if (n === 0) {
 return sign * (0, _j.y0)(x);
 }
 if (n === 1) {
 return sign * (0, _j2.y1)(x);
 }

 /* test for overflow */
 if (x === 0.0) {
 // mtherr('yn', SING);
 return -Infinity * sign;
 } else if (x < 0.0) {
 // mtherr('yn', DOMAIN);
 return NaN;
 }

 /* forward recurrence on n */
 anm2 = (0, _j.y0)(x);
 anm1 = (0, _j2.y1)(x);
 k = 1;
 r = 2 * k;
 do {
 an = r * anm1 / x - anm2;
 anm2 = anm1;
 anm1 = an;
 r += 2.0;
 ++k;
 } while (k < n);

 return sign * an;
}

exports.yn = yn;
},{"./j0.js":67,"./j1.js":68}],81:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.yv = undefined;

var _jv = require('./jv.js');

var _yn = require('./yn.js');

/**
 * @file yv.js Bessel function of second kind of non-integer order
 *
 * Cephes Math Library Release 2.8: June, 2000
 * Copyright 1984, 1987, 2000 by Stephen L. Moshier
 * Ported to ECMAScript 2018
 * Copyright (c) 2018, Kings Distributed Systems
 *
 * @author KC Erb, kc@kcerb.com
 * @date April 2018
 */
function yv(v, x) {
 var y = void 0,
 t = void 0,
 n = void 0;
 n = Math.trunc(v);
 if (n === v) {
 y = (0, _yn.yn)(n, x);
 return y;
 } else if (v === Math.floor(v)) {
 /* Zero in denominator. */
 // mtherr('yv', DOMAIN);
 return NaN;
 }

 t = Math.PI * v;
 y = (Math.cos(t) * (0, _jv.jv)(v, x) - (0, _jv.jv)(-v, x)) / Math.sin(t);

 if (!isFinite(y)) {
 if (v > 0) {
 // mtherr('yv', OVERFLOW);
 return -Infinity;
 } else if (v < -1e10) {
 /* Whether it's +inf or -inf is numerically ill-defined. */
 // mtherr('yv', DOMAIN);
 return NaN;
 }
 }

 return y;
}

exports.yv = yv;
},{"./jv.js":69,"./yn.js":80}],82:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});

var _createClass = function () { function defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } } return function (Constructor, protoProps, staticProps) { if (protoProps) defineProperties(Constructor.prototype, protoProps); if (staticProps) defineProperties(Constructor, staticProps); return Constructor; }; }();

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

var Chebyshev = function () {
 function Chebyshev() {
 _classCallCheck(this, Chebyshev);
 }

 _createClass(Chebyshev, null, [{
 key: 'tn',

 // SOURCE: https://en.wikipedia.org/wiki/Chebyshev_polynomials#Explicit_expressions
 // TODO: the following only accepts integer `n`s -
 // Make continuous via the Gauss hypergeometric function
 // https://github.com/scipy/scipy/blob/master/scipy/special/cephes/hyp2f1.c
 value: function tn(n, viewX, x, viewR, r) {
 var arg = void 0;
 switch (true) {
 case Math.abs(viewX[x]) <= 1:
 arg = n * Math.acos(viewX[x]);
 viewR[r] = Math.cos(arg);
 break;
 case viewX[x] > 1:
 arg = n * Math.acosh(viewX[x]);
 viewR[r] = Math.cosh(arg);
 break;
 default:
 var sign = Math.pow(-1, n);
 arg = n * Math.acosh(-viewX[x]);
 viewR[r] = sign * Math.cosh(arg);
 }
 }

 // SOURCE: https://en.wikipedia.org/wiki/Chebyshev_polynomials#Explicit_expressions
 // TODO: same as tn

 }, {
 key: 'un',
 value: function un(n, viewX, x, viewR, r) {
 switch (true) {
 case Math.abs(viewX[x]) > 1:
 var rad = Math.sqrt(viewX[x] * viewX[x] - 1);
 var numerator = Math.pow(viewX[x] + rad, n + 1) - Math.pow(viewX[x] - rad, n + 1);
 var denominator = 2 * rad;
 viewR[r] = numerator / denominator;
 break;
 case Math.abs(viewX[x]) <= 1:
 // TODO, implement as hypergeometric to cover all ranges
 break;
 default:
 }
 }

 // Approximate the function `func` in the interval [a, b] with
 // an order-n Chebyshev polynomial (T).
 // f ≈ ∑c_k T_k(y) - c_0/2 where k = 0->(m-1) in integer steps.
 // Return the coeffecients c_n of that polynomial.
 // TODO - Numerical Recipes suggests that if `Math.cos` is bogging this down,
 // we should look at 12.3 and consider "fast cosine transform methods"

 }, {
 key: 'approx',
 value: function approx(func, a, b) {
 var n = arguments.length > 3 && arguments[3] !== undefined ? arguments[3] : 50;

 var y = void 0;
 var c = new Float64Array(n);
 var f = new Float64Array(n);

 var bma = 0.5 * (b - a);
 var bpa = 0.5 * (b + a);

 // Evaluate func at n points
 for (var k = 0; k < n; k++) {
 y = Math.cos(Math.PI * (k + 0.5) / n);
 f[k] = func(y * bma + bpa);
 }

 // Calculate the c_j
 var fac = 2.0 / n;
 var sum = 0.0;
 for (var j = 0; j < n; j++) {
 sum = 0.0;
 for (var _k = 0; _k < n; _k++) {
 sum += f[_k] * Math.cos(Math.PI * j * (_k + 0.5) / n);
 }

 c[j] = fac * sum;
 }

 return c;
 }

 // Chebyshev evaluation: The Chebyshev polynomial ∑c_k T_k(y) - c_0/2 where k=0->(m-1)
 // is evaluated at y = [x - (b+a)/2] / [(b-a)/2] using the coeffecients c. These coeffecients
 // can be calculated by AdvMathChebyshev::approx.

 }, {
 key: 'eval',
 value: function _eval(viewX, x, viewR, r, a, b, c) {
 var opts = arguments.length > 7 && arguments[7] !== undefined ? arguments[7] : {};

 var outOfRange = (viewX[x] - a) * (viewX[x] - b) > 0.0;
 if (outOfRange) {
 throw new Error('x not in range in Chebyshev::eval');
 }

 var sv = void 0,
 j = void 0;
 var d = 0.0;
 var dd = 0.0;

 // Change of variable
 var y = (2.0 * viewX[x] - a - b) / (b - a);
 var y2 = 2.0 * y;

 // Clenshaw’s recurrence
 var m = this.setm(c, opts);
 for (j = m - 1; j > 0; j--) {
 sv = d;
 d = y2 * d - dd + c[j];
 dd = sv;
 }

 viewR[r] = y * d - dd + 0.5 * c[0];
 }

 // Don't bother iterating over coeffecients less than opts.thresh

 }, {
 key: 'setm',
 value: function setm(c, opts) {
 var m = c.length;

 if (opts.thresh) {
 while (m > 1 && Math.abs(c[m - 1]) < opts.thresh) {
 m--;
 }
 } else if (opts.m) {
 m = opts.m;
 }
 return m;
 }
 }]);

 return Chebyshev;
}();

exports.default = Chebyshev;
},{}],83:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.sinpi = exports.rgamma = exports.gamma = undefined;

var _slicedToArray = function () { function sliceIterator(arr, i) { var _arr = []; var _n = true; var _d = false; var _e = undefined; try { for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) { _arr.push(_s.value); if (i && _arr.length === i) break; } } catch (err) { _d = true; _e = err; } finally { try { if (!_n && _i["return"]) _i["return"](); } finally { if (_d) throw _e; } } return _arr; } return function (arr, i) { if (Array.isArray(arr)) { return arr; } else if (Symbol.iterator in Object(arr)) { return sliceIterator(arr, i); } else { throw new TypeError("Invalid attempt to destructure non-iterable instance"); } }; }();

var _complex = require('../../utils/complex.js');

var _complex2 = _interopRequireDefault(_complex);

var _pythonHelpers = require('../../utils/pythonHelpers.js');

var py = _interopRequireWildcard(_pythonHelpers);

function _interopRequireWildcard(obj) { if (obj && obj.__esModule) { return obj; } else { var newObj = {}; if (obj != null) { for (var key in obj) { if (Object.prototype.hasOwnProperty.call(obj, key)) newObj[key] = obj[key]; } } newObj.default = obj; return newObj; } }

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

var EXACT_GAMMA = [Infinity, 1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320.0, 362880.0, 3628800.0, 39916800.0, 479001600.0, 6227020800.0, 87178291200.0, 1307674368000.0, 20922789888000.0, 355687428096000.0, 6402373705728000.0, 121645100408832000.0, 2432902008176640000.0];

// Lanczos coefficients used by the GNU Scientific Library
var LANCZOS_G = 7;
var LANCZOS_P = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];

function gammaReal(x) {
 var _intx = Math.round(x);
 if (_intx === x) {
 if (_intx <= 0) {
 // return (-1)**_intx * INF
 throw new ZeroDivisionError('gamma function pole');
 }
 if (EXACT_GAMMA[_intx]) {
 return EXACT_GAMMA[_intx];
 }
 }

 if (x < 0.5) {
 // TODO: sinpi
 return Math.PI / (sinpiReal(x) * gammaReal(1 - x));
 } else {
 x -= 1.0;
 var r = LANCZOS_P[0];
 for (var i = 1; i < LANCZOS_G + 2; i++) {
 r += LANCZOS_P[i] / (x + i);
 }
 var t = x + LANCZOS_G + 0.5;
 return 2.506628274631000502417 * t ** (x + 0.5) * Math.exp(-t) * r;
 }
}

function gammaComplex(x) {
 if (_complex2.default.re(x) < 0.5) {
 // TODO: sinpi
 var oneMinusX = _complex2.default.sub(1, x);
 var denom = _complex2.default.mul(sinpiComplex(x), gammaComplex(oneMinusX));
 return _complex2.default.div(Math.PI, denom);
 } else {
 x = _complex2.default.sub(x, 1.0);
 var r = LANCZOS_P[0];
 for (var i = 1; i < LANCZOS_G + 2; i++) {
 var xPlusI = _complex2.default.add(x, i);
 var term = _complex2.default.div(LANCZOS_P[i], xPlusI);
 r = _complex2.default.add(r, term);
 }
 var t = _complex2.default.add(x, LANCZOS_G + 0.5);
 // return 2.506628274631000502417 * t**(x+0.5) * Complex.exp(-t) * r
 var xPlusHalf = _complex2.default.add(x, 0.5);
 var tToTheXPlusHalf = _complex2.default.pow(t, xPlusHalf);
 var negT = _complex2.default.mul(-1, t);
 var rTimesEToTheMinusT = _complex2.default.mul(r, _complex2.default.exp(negT));
 var product = _complex2.default.mul(tToTheXPlusHalf, rTimesEToTheMinusT);
 return _complex2.default.mul(2.506628274631000502417, product);
 }
}

function gamma(x) {
 if (x.constructor === Array) {
 return gammaComplex(x);
 } else {
 return gammaReal(x);
 }
}

function rgamma(x) {
 try {
 return _complex2.default.inverse(gamma(x));
 } catch (e) {
 if (e.name === 'ZeroDivisionError') {
 return _complex2.default.mul(x, 0.0);
 } else {
 // only catch ZeroDivisionError
 throw e;
 }
 }
}

function sinpiReal(x) {
 if (x < 0) return -sinpiReal(-x);

 var _py$divmod = py.divmod(x, 0.5),
 _py$divmod2 = _slicedToArray(_py$divmod, 2),
 n = _py$divmod2[0],
 r = _py$divmod2[1];

 r *= Math.PI;
 n = py.mod(n, 4);
 if (n === 0) return Math.sin(r);
 if (n === 1) return Math.cos(r);
 if (n === 2) return -Math.sin(r);
 if (n === 3) return -Math.cos(r);
}

function sinpiComplex(z) {
 if (_complex2.default.re(z) < 0) {
 var negZ = _complex2.default.mul(-1, z);
 var res = sinpiComplex(negZ);
 return _complex2.default.mul(-1, res);
 };

 var _py$divmod3 = py.divmod(_complex2.default.re(z), 0.5),
 _py$divmod4 = _slicedToArray(_py$divmod3, 2),
 n = _py$divmod4[0],
 r = _py$divmod4[1];

 z = _complex2.default.mul(Math.PI, [r, _complex2.default.im(z)]);
 n = py.mod(n, 4);
 if (n === 0) return _complex2.default.sin(z);
 if (n === 1) return _complex2.default.cos(z);
 if (n === 2) return _complex2.default.mul(-1, _complex2.default.sin(z));
 if (n === 3) return _complex2.default.mul(-1, _complex2.default.cos(z));
}

function sinpi(x) {
 if (x.constructor === Array) {
 return sinpiComplex(x);
 } else {
 return sinpiReal(x);
 }
}

function ZeroDivisionError(message) {
 this.message = message;
 this.name = 'ZeroDivisionError';
}

exports.gamma = gamma;
exports.rgamma = rgamma;
exports.sinpi = sinpi;
},{"../../utils/complex.js":89,"../../utils/pythonHelpers.js":94}],84:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.reset = exports.j = exports.expj = exports.expjpi = exports._defaultHyperMaxprec = exports.eps = exports._fixedPrecision = exports.dps = exports.prec = exports.isint = exports.NoConvergence = exports.isComplexType = exports._isRealType = exports.nstr = exports.convertParam = exports.hypsum = exports.fprod = exports.isnpint = exports.mag = exports.nintDistance = exports.convert = undefined;

var _pythonHelpers = require('../../utils/pythonHelpers.js');

var py = _interopRequireWildcard(_pythonHelpers);

var _complex = require('../../utils/complex.js');

var _complex2 = _interopRequireDefault(_complex);

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _interopRequireWildcard(obj) { if (obj && obj.__esModule) { return obj; } else { var newObj = {}; if (obj != null) { for (var key in obj) { if (Object.prototype.hasOwnProperty.call(obj, key)) newObj[key] = obj[key]; } } newObj.default = obj; return newObj; } }

// Make fixed precision
var prec = 53;
var dps = 15;
var _fixedPrecision = true;
var eps = Number.EPSILON;
var j = [0, 1];

function reset() {
 exports.prec = prec = 53;
}

// mpmath relies on a bunch of helper functions stored in `ctx`.
function nintDistance(z) {
 var n = void 0,
 distance = void 0;
 if (z.constructor === Array) {
 n = Math.round(z[0]);
 } else {
 n = Math.round(z);
 }

 if (n === z) {
 distance = -Infinity;
 } else {
 distance = mag(_complex2.default.sub(z, n));
 }

 return [n, distance];
}

// function magc(z) {
// return Math.max(mag(Complex.re(z)), mag(Complex.im(z))) + 1;
// }

function mag(x) {
 if (x) {
 return py.frexp(_complex2.default.abs(x))[1];
 } else {
 return -Infinity;
 }
}

function isnpint(x) {
 if (x.constructor === Array) {
 if (x[1] !== 0) {
 return false;
 }
 x = x[0];
 }
 return x <= 0.0 && Math.round(x) === x;
}

function fprod(args) {
 var prod = 1;
 var _iteratorNormalCompletion = true;
 var _didIteratorError = false;
 var _iteratorError = undefined;

 try {
 for (var _iterator = args[Symbol.iterator](), _step; !(_iteratorNormalCompletion = (_step = _iterator.next()).done); _iteratorNormalCompletion = true) {
 var arg = _step.value;

 prod = _complex2.default.mul(prod, arg);
 }
 } catch (err) {
 _didIteratorError = true;
 _iteratorError = err;
 } finally {
 try {
 if (!_iteratorNormalCompletion && _iterator.return) {
 _iterator.return();
 }
 } finally {
 if (_didIteratorError) {
 throw _iteratorError;
 }
 }
 }

 return prod;
}

// function hypsumMp(p, q, flags, coeffs, z, opts={}){
// if (!('accurate_small' in opts)) {
// opts['accurate_small'] = true;
// }
//
// if (hasattr(z, "_mpf_")) {
// key = p, q, flags, 'R';
// v = z._mpf_;
// } else if (hasattr(z, "_mpc_")) {
// key = p, q, flags, 'C';
// v = z._mpc_;
// }
// if (key not in ctx.hyp_summators) {
// ctx.hyp_summators[key] = libmp.make_hyp_summator(key)[1];
// };
// summator = ctx.hyp_summators[key];
// prec = ctx.prec;
// maxprec = kwargs.get('maxprec', ctx._default_hyper_maxprec(prec));
// extraprec = 50;
// epsshift = 25;
// // Jumps in magnitude occur when parameters are close to negative
// // integers. We must ensure that these terms are included in
// // the sum and added accurately
// magnitude_check = {};
// max_total_jump = 0;
// for (i, c in enumerate(coeffs)) {
// if (flags[i]==='Z') {
// if (i >= p and c <= 0) {
// ok = False;
// for (ii, cc in enumerate(coeffs[:p])) {
// // Note: c <= cc or c < cc, depending on convention
// if (flags[ii]==='Z' and cc <= 0 and c <= cc) {
// ok = True;
// }
// }
// if (not ok) {
// raise ZeroDivisionError("pole in hypergeometric series")
// }
// }
// continue
// }
// [n, d] = ctx.nint_distance(c);
// n = -int(n);
// d = -d;
// if (i >= p and n >= 0 and d > 4) {
// if (n in magnitude_check) {
// magnitude_check[n] += d;
// } else {
// magnitude_check[n] = d
// }
// extraprec = max(extraprec, d - prec + 60)
// }
// max_total_jump += abs(d);
// }
// while (1) {
// if (extraprec > maxprec) {
// raise ValueError(ctx._hypsum_msg % (prec, prec+extraprec))
// }
// wp = prec + extraprec;
// if (magnitude_check) {
// mag_dict = dict((n,None) for n in magnitude_check);
// } else {
// mag_dict = {}
// }
// [zv, have_complex, magnitude] = summator(coeffs, v, prec, wp, epsshift, mag_dict, **kwargs);
// cancel = -magnitude;
// jumps_resolved = true;
// if (extraprec < max_total_jump) {
// for (n in mag_dict.values()) {
// if ((n is None) or (n < prec)) {
// jumps_resolved = false;
// break
// }
// }
// }
// accurate = (cancel < extraprec-25-5 or not accurate_small)
// if (jumps_resolved) {
// if (accurate) {
// break
// }
// //zero?
// zeroprec = kwargs.get('zeroprec');
// if (zeroprec is not None) {
// if (cancel > zeroprec) {
// if (have_complex) {
// return ctx.mpc(0);
// } else {
// return ctx.zerol
// }
// }
// }
// }
// // Some near-singularities were not included, so increase
// // precision and repeat until they are
// extraprec *= 2;
// // Possible workaround for bad roundoff in fixed-point arithmetic
// epsshift += 5;
// extraprec += 5;
// }
// if (type(zv) is tuple) {
// if (have_complex) {
// return ctx.make_mpc(zv);
// } else {
// return ctx.make_mpf(zv);
// }
// } else {
// return zv;
// }
// }

// Depending on the context, mpmath might call any one of a number of hypsum functions
// this one got translated by accident for our immediate usecase, but may be useful later
function hypsum(p, q, types, coeffs, z) {
 var opts = arguments.length > 5 && arguments[5] !== undefined ? arguments[5] : {};

 if (!('maxterms' in opts)) {
 opts['maxterms'] = 6000;
 }

 var num = py.range(p);
 var den = py.range(p, p + q);
 var tol = eps;
 var s = 1.0;
 var t = 1.0;
 var k = 0;

 while (1) {
 var _iteratorNormalCompletion2 = true;
 var _didIteratorError2 = false;
 var _iteratorError2 = undefined;

 try {
 for (var _iterator2 = num[Symbol.iterator](), _step2; !(_iteratorNormalCompletion2 = (_step2 = _iterator2.next()).done); _iteratorNormalCompletion2 = true) {
 var i = _step2.value;
 t = _complex2.default.mul(t, _complex2.default.add(coeffs[i], k));
 }
 } catch (err) {
 _didIteratorError2 = true;
 _iteratorError2 = err;
 } finally {
 try {
 if (!_iteratorNormalCompletion2 && _iterator2.return) {
 _iterator2.return();
 }
 } finally {
 if (_didIteratorError2) {
 throw _iteratorError2;
 }
 }
 }

 var _iteratorNormalCompletion3 = true;
 var _didIteratorError3 = false;
 var _iteratorError3 = undefined;

 try {
 for (var _iterator3 = den[Symbol.iterator](), _step3; !(_iteratorNormalCompletion3 = (_step3 = _iterator3.next()).done); _iteratorNormalCompletion3 = true) {
 var _i = _step3.value;
 t = _complex2.default.div(t, _complex2.default.add(coeffs[_i], k));
 }
 } catch (err) {
 _didIteratorError3 = true;
 _iteratorError3 = err;
 } finally {
 try {
 if (!_iteratorNormalCompletion3 && _iterator3.return) {
 _iterator3.return();
 }
 } finally {
 if (_didIteratorError3) {
 throw _iteratorError3;
 }
 }
 }

 k += 1;
 t = _complex2.default.div(t, k);
 t = _complex2.default.mul(t, z);
 s = _complex2.default.add(s, t);
 if (_complex2.default.abs(t) < tol) {
 return s;
 }
 if (k > opts.maxterms) {
 throw new NoConvergence('hypsum failed to converge');
 }
 }
}

// original divided when given a tuple, since js doesn't have tuples or
// much in the way of numeric types, we'll leave that out and assume numbers
// passed are either ints, real, or complex.
function convertParam(z) {
 // if (tuple) {
 // [p, q] = z;
 // return [p / q, 'R'];
 // }
 var intz = void 0;
 if (z.constructor === Array) {
 intz = Math.floor(z[0]);
 } else {
 intz = Math.floor(z);
 }

 if (z === intz) {
 return [intz, 'Z'];
 }

 return [z, 'R'];
}

function convert(x) {
 return _complex2.default.ensureComplex(x);
}

function nstr(x) {
 var n = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : 5;
 var opts = arguments.length > 2 && arguments[2] !== undefined ? arguments[2] : {};

 // x = Complex.ensureComplex(x);

 // TODO - to string function - mpi
 // if (hasattr(x, "_mpi_")) {
 // return libmp.mpi_to_str(x._mpi_, n, opts)
 // }
 // if (hasattr(x, "_mpci_")) {
 // re = libmp.mpi_to_str(x._mpci_[0], n, opts)
 // im = libmp.mpi_to_str(x._mpci_[1], n, opts)
 // return "(%s + %s*j)" % (re, im)
 // }
 return String(x);
}

function isint(z) {
 if (z.constructor === Array) {
 if (_complex2.default.im(z)) {
 return false;
 }
 z = _complex2.default.re(z);
 }
 return Number.isInteger(z);
}

function expjpi(x) {
 return _complex2.default.exp(_complex2.default.mul([0, Math.PI], x));
}

function expj(x) {
 return _complex2.default.exp(_complex2.default.mul([0, 1], x));
}

function _isRealType(z) {
 return z.constructor !== Array;
}

function isComplexType(z) {
 return z.constructor === Array;
}

function NoConvergence(message) {
 this.message = message;
 this.name = 'NoConvergence';
}

function _defaultHyperMaxprec(p) {
 return Math.trunc(1000 * p ** 0.25 + 4 * p);
}

exports.convert = convert;
exports.nintDistance = nintDistance;
exports.mag = mag;
exports.isnpint = isnpint;
exports.fprod = fprod;
exports.hypsum = hypsum;
exports.convertParam = convertParam;
exports.nstr = nstr;
exports._isRealType = _isRealType;
exports.isComplexType = isComplexType;
exports.NoConvergence = NoConvergence;
exports.isint = isint;
exports.prec = prec;
exports.dps = dps;
exports._fixedPrecision = _fixedPrecision;
exports.eps = eps;
exports._defaultHyperMaxprec = _defaultHyperMaxprec;
exports.expjpi = expjpi;
exports.expj = expj;
exports.j = j;
exports.reset = reset;
},{"../../utils/complex.js":89,"../../utils/pythonHelpers.js":94}],85:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.struvel = exports.besselk = exports.besseli = undefined;

var _complex = require('../../../utils/complex.js');

var _complex2 = _interopRequireDefault(_complex);

var _hypergeometric = require('./hypergeometric.js');

var mpHypergeometric = _interopRequireWildcard(_hypergeometric);

var _ctx = require('../ctx.js');

var ctx = _interopRequireWildcard(_ctx);

function _interopRequireWildcard(obj) { if (obj && obj.__esModule) { return obj; } else { var newObj = {}; if (obj != null) { for (var key in obj) { if (Object.prototype.hasOwnProperty.call(obj, key)) newObj[key] = obj[key]; } } newObj.default = obj; return newObj; } }

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

// besselI of complex order and arg of selected derivative
function besseli(n, z) {
 var derivative = arguments.length > 2 && arguments[2] !== undefined ? arguments[2] : 0;
 var opts = arguments.length > 3 && arguments[3] !== undefined ? arguments[3] : {};

 ctx.reset();
 var v = void 0,
 r = void 0,
 t = void 0;
 n = ctx.convert(n);
 z = ctx.convert(z);

 if (_complex2.default.isZero(z)) {
 if (derivative !== 0) {
 // raise ValueError
 throw new Error('besseli derivative not defined for z[0] = z[1] = 0');
 }

 if (_complex2.default.isZero(n)) {
 // I(0,0) = 1
 return 1;
 }

 // Integer orders are zero when z is zero
 if (n[1] === 0 && Number.isInteger(n[0])) {
 return 0;
 }

 // Non-integer orders (the real part) are either NaN, 0, or Infinity
 if (n[0] === 0) {
 return NaN;
 } else if (n[0] > 0) {
 return 0;
 } else {
 return Infinity;
 }
 }

 if (derivative) {
 var h = function h(n, d) {
 var zSquared = _complex2.default.mul(z, z);
 r = _complex2.default.mul(zSquared, 0.25);
 var nOver2 = _complex2.default.div(n, 2);
 var nPlus1 = _complex2.default.add(n, 1);
 var b = [_complex2.default.add(nOver2, 0.5 - d / 2), _complex2.default.add(nOver2, 1 - d / 2), nPlus1];
 // mpmath returns a tuple wrapped in a list
 // [([2,ctx.pi,z],[d-2*n,0.5,n-d],[n+1],B,[(n+1)*0.5,(n+2)*0.5],B,r)]
 t = [[2, Math.PI, z], // ws
 // [d-2*n,0.5,n-d]
 [_complex2.default.sub(d, _complex2.default.mul(2, n)), 0.5, _complex2.default.sub(n, d)], // cs
 [nPlus1], // alphas
 b, // betas
 // [(n+1)*0.5,(n+2)*0.5]
 [_complex2.default.mul(nPlus1, 0.5), _complex2.default.mul(_complex2.default.add(n, 2), 0.5)], // as = 2
 b, // bs = 1
 r // z
 ];
 return [t];
 };
 v = mpHypergeometric.hypercomb(h, [n, derivative], opts);
 } else {
 var _h = function _h(n) {
 var nPlus1 = _complex2.default.add(n, 1);
 var w = _complex2.default.mul(z, 0.5);
 r = _complex2.default.mul(w, w);
 // as = 0 bs = 1
 t = [[w], [n], [], [nPlus1], [], [nPlus1], r];
 return [t];
 };
 v = mpHypergeometric.hypercomb(_h, [n], opts);
 }
 return v;
}

// complex arg - besselK
function besselk(n, z) {
 var opts = arguments.length > 2 && arguments[2] !== undefined ? arguments[2] : {};

 ctx.reset();
 n = _complex2.default.ensureComplex(n);
 z = _complex2.default.ensureComplex(z);
 if (_complex2.default.isZero(z)) {
 return Infinity;
 }

 var h = void 0;
 var capitalM = ctx.mag(z);
 if (capitalM < 1) {
 // Represent as limit definition
 h = function h(n) {
 var nPlus1 = _complex2.default.add(n, 1);
 var zOver2 = _complex2.default.div(z, 2);
 var negN = _complex2.default.mul(n, -1);
 var r = _complex2.default.pow(zOver2, 2);
 // ws cs alphas betas as=0 bs=1
 var t1 = [[z, 2], [negN, _complex2.default.sub(n, 1)], [n], [], [], [_complex2.default.sub(1, n)], r];
 var t2 = [[z, 2], [n, _complex2.default.mul(nPlus1, -1)], [negN], [], [], [nPlus1], r];
 return [t1, t2];
 };
 // We could use the limit definition always, but it leads
 // to very bad cancellation (of exponentially large terms)
 // for large real z
 // Instead represent in terms of 2F0
 } else {
 ctx.prec += capitalM;
 h = function h(n) {
 var negZ = _complex2.default.mul(z, -1);
 var t = [[Math.PI / 2, z, _complex2.default.exp(negZ)], [0.5, -0.5, 1], [], [], [_complex2.default.add(n, 0.5), _complex2.default.sub(0.5, n)], // a=2
 [], // b=0
 _complex2.default.mul(-0.5, _complex2.default.inverse(z))];
 return [t];
 };
 }

 return mpHypergeometric.hypercomb(h, [n], opts);
}

function struvel(n, z, opts) {
 ctx.reset();
 n = ctx.convert(n);
 z = ctx.convert(z);
 // http://functions.wolfram.com/Bessel-TypeFunctions/StruveL/26/01/02/
 var h = function h(n) {
 var capitalT = [[_complex2.default.div(z, 2), 0.5 * Math.sqrt(Math.PI)], [_complex2.default.add(n, 1), -1], [], [_complex2.default.add(n, 1.5)], [1], [1.5, _complex2.default.add(n, 1.5)], _complex2.default.pow(_complex2.default.div(z, 2), 2)];
 return [capitalT];
 };
 return mpHypergeometric.hypercomb(h, [n], opts);
}

exports.besseli = besseli;
exports.besselk = besselk;
exports.struvel = struvel;
},{"../../../utils/complex.js":89,"../ctx.js":84,"./hypergeometric.js":88}],86:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.gammaprod = undefined;

var _ctx = require('../ctx.js');

var ctx = _interopRequireWildcard(_ctx);

var _complex = require('../../../utils/complex.js');

var _complex2 = _interopRequireDefault(_complex);

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _interopRequireWildcard(obj) { if (obj && obj.__esModule) { return obj; } else { var newObj = {}; if (obj != null) { for (var key in obj) { if (Object.prototype.hasOwnProperty.call(obj, key)) newObj[key] = obj[key]; } } newObj.default = obj; return newObj; } }

function gammaprod(a, b) {
 var _infsign = arguments.length > 2 && arguments[2] !== undefined ? arguments[2] : false;

 a = a.map(function (x) {
 return ctx.convert(x);
 });
 b = b.map(function (x) {
 return ctx.convert(x);
 });
 var polesNum = [];
 var polesDen = [];
 var regularNum = [];
 var regularDen = [];

 var sign = void 0;

 a.map(function (x) {
 if (ctx.isnpint(x)) {
 polesNum.push(x);
 } else {
 regularNum.push(x);
 }
 });

 b.map(function (x) {
 if (ctx.isnpint(x)) {
 polesDen.push(x);
 } else {
 regularDen.push(x);
 }
 });

 // One more pole in numerator or denominator gives 0 or inf
 if (polesNum.length < polesDen.length) return 0;
 if (polesNum.length > polesDen.length) {
 // Get correct sign of infinity for x+h, h -> 0 from above
 // XXX: hack, this should be done properly
 if (_infsign) {
 var _a, _b;

 a = polesNum.map(function (x) {
 return x && (_complex2.default.mul(x, 1 + ctx.eps) || _complex2.default.add(x, ctx.eps));
 });
 b = polesDen.map(function (x) {
 return x && (_complex2.default.mul(x, 1 + ctx.eps) || _complex2.default.add(x, ctx.eps));
 });
 sign = ctx.sign(ctx.gammaprod((_a = a).push.apply(_a, regularNum), (_b = b).push.apply(_b, regularDen)));
 return _complex2.default.mul(sign, Infinity);
 } else {
 return Infinity;
 }
 }
 // All poles cancel
 // lim G(i)/G(j) = (-1)**(i+j) * gamma(1-j) / gamma(1-i)
 var p = 1;
 var orig = ctx.prec;
 try {
 ctx.prec = orig + 15;
 while (polesNum.length > 0) {
 var i = polesNum.pop();
 var j = polesDen.pop();
 var iPlusJ = _complex2.default.add(i, j);
 sign = _complex2.default.pow(-1, iPlusJ);
 var gammaOneMinusJ = ctx.gamma(_complex2.default.sub(1, j));
 var gammaOneMinusI = ctx.gamma(_complex2.default.sub(1, i));
 var quotient = _complex2.default.div(gammaOneMinusJ, gammaOneMinusI);
 p *= _complex2.default.mul(sign, quotient);
 }
 regularNum.forEach(function (x) {
 p = _complex2.default.mul(p, ctx.gamma(x));
 });
 regularDen.forEach(function (x) {
 p = _complex2.default.div(p, ctx.gamma(x));
 });
 } finally {
 ctx.prec = orig;
 }
 return p;
}

exports.gammaprod = gammaprod;
},{"../../../utils/complex.js":89,"../ctx.js":84}],87:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.sign = undefined;

var _complex = require('../../../utils/complex.js');

var _complex2 = _interopRequireDefault(_complex);

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function sign(ctx, x) {
 if (!x || Number.isNaN(x)) return x;
 if (ctx._isRealType(x)) {
 if (x > 0) {
 return 1;
 } else {
 return -1;
 }
 }

 return _complex2.default.div(x, _complex2.default.abs(x));
}

exports.sign = sign;
},{"../../../utils/complex.js":89}],88:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.hyp1f2 = exports.hyp2f0 = exports.hyper = exports.hypercomb = undefined;

var _slicedToArray = function () { function sliceIterator(arr, i) { var _arr = []; var _n = true; var _d = false; var _e = undefined; try { for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) { _arr.push(_s.value); if (i && _arr.length === i) break; } } catch (err) { _d = true; _e = err; } finally { try { if (!_n && _i["return"]) _i["return"](); } finally { if (_d) throw _e; } } return _arr; } return function (arr, i) { if (Array.isArray(arr)) { return arr; } else if (Symbol.iterator in Object(arr)) { return sliceIterator(arr, i); } else { throw new TypeError("Invalid attempt to destructure non-iterable instance"); } }; }();

var _complex = require('../../../utils/complex.js');

var _complex2 = _interopRequireDefault(_complex);

var _ctx = require('../ctx.js');

var ctx = _interopRequireWildcard(_ctx);

var _pythonHelpers = require('../../../utils/pythonHelpers.js');

var py = _interopRequireWildcard(_pythonHelpers);

var _math = require('../../math2/math2.js');

var math2 = _interopRequireWildcard(_math);

var _factorials = require('./factorials.js');

var mpFactorials = _interopRequireWildcard(_factorials);

var _functions = require('./functions.js');

var mpFunctions = _interopRequireWildcard(_functions);

function _interopRequireWildcard(obj) { if (obj && obj.__esModule) { return obj; } else { var newObj = {}; if (obj != null) { for (var key in obj) { if (Object.prototype.hasOwnProperty.call(obj, key)) newObj[key] = obj[key]; } } newObj.default = obj; return newObj; } }

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _toConsumableArray(arr) { if (Array.isArray(arr)) { for (var i = 0, arr2 = Array(arr.length); i < arr.length; i++) { arr2[i] = arr[i]; } return arr2; } else { return Array.from(arr); } }

// In original this is called hypercomb

// Thoughts:
// I'd rather it was given / returned a keyed object.

function hypercomb(func) {
 var params = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : [];
 var opts = arguments.length > 2 && arguments[2] !== undefined ? arguments[2] : {};

 var hmag = void 0,
 zeroOk = void 0,
 infOk = void 0;
 var orig = ctx.prec;
 var discardKnownZeros = py.get(opts, 'discardKnownZeros', true);
 var sumvalue = 0;
 var origParams = [].concat(_toConsumableArray(params));
 var verbose = py.get(opts, 'verbose', false);
 var maxprec = py.get(opts, 'maxprec', ctx._defaultHyperMaxprec(orig));
 opts.maxprec = maxprec; // For calls to hypsum
 var zeroprec = py.get(opts, 'zeroprec');
 var infprec = py.get(opts, 'infprec');
 var perturbedReferenceValue = null;
 var hextra = 0;
 try {
 while (1) {
 ctx.prec += 10;
 if (ctx.prec > maxprec) {
 // raise ValueError(_hypercomb_msg % (orig, ctx.prec));
 throw new Error('_hypercomb_msg: orig: ' + orig + ' | ctx.prec: ' + ctx.prec);
 }

 var orig2 = ctx.prec;
 params = [].concat(_toConsumableArray(origParams));
 var terms = func.apply(undefined, _toConsumableArray(params));
 if (verbose) {
  console.log();
  console.log('ENTERING hypercomb main loop');
  console.log('prec =', ctx.prec);
  console.log('hextra', hextra);
 }

 var _checkNeedPerturb2 = _checkNeedPerturb(terms, orig, discardKnownZeros),
 _checkNeedPerturb3 = _slicedToArray(_checkNeedPerturb2, 4),
 perturb = _checkNeedPerturb3[0],
 recompute = _checkNeedPerturb3[1],
 extraprec = _checkNeedPerturb3[2],
 discard = _checkNeedPerturb3[3];

 ctx.prec += extraprec;
 if (perturb) {
 if (opts.hasOwnProperty('hmag')) {
 hmag = opts.hmag;
 } else if (ctx._fixedPrecision) {
 hmag = Math.trunc(ctx.prec * 0.3);
 } else {
 hmag = orig + 10 + hextra;
 }
 opts.perturbParams = [];
 var h = py.ldexp(1, -hmag);
 ctx.prec = orig2 + 10 + hmag + 10;
 var _iteratorNormalCompletion = true;
 var _didIteratorError = false;
 var _iteratorError = undefined;

 try {
 for (var _iterator = py.range(params.length)[Symbol.iterator](), _step; !(_iteratorNormalCompletion = (_step = _iterator.next()).done); _iteratorNormalCompletion = true) {
 var k = _step.value;

 opts.perturbParams[k] = h;
 params[k] = _complex2.default.add(params[k], h);
 // Heuristically ensure that the perturbations
 // are "independent" so that two perturbations
 // don't accidentally cancel each other out
 // in a subtraction.
 h += h / (k + 1);
 }
 } catch (err) {
 _didIteratorError = true;
 _iteratorError = err;
 } finally {
 try {
 if (!_iteratorNormalCompletion && _iterator.return) {
 _iterator.return();
 }
 } finally {
 if (_didIteratorError) {
 throw _iteratorError;
 }
 }
 }
 }
 if (recompute) {
 terms = func.apply(undefined, _toConsumableArray(params).concat([opts]));
 }
 if (discardKnownZeros) {
 var _iteratorNormalCompletion2 = true;
 var _didIteratorError2 = false;
 var _iteratorError2 = undefined;

 try {
 for (var _iterator2 = discard[Symbol.iterator](), _step2; !(_iteratorNormalCompletion2 = (_step2 = _iterator2.next()).done); _iteratorNormalCompletion2 = true) {
 var index = _step2.value;
 delete terms[index];
 }
 } catch (err) {
 _didIteratorError2 = true;
 _iteratorError2 = err;
 } finally {
 try {
 if (!_iteratorNormalCompletion2 && _iterator2.return) {
 _iterator2.return();
 }
 } finally {
 if (_didIteratorError2) {
 throw _iteratorError2;
 }
 }
 }

 terms = terms.filter(Boolean);
 }
 if (!terms.length) {
 return 0;
 }

 // Now actually perform the calculation!
 var evaluatedTerms = [];
 var _iteratorNormalCompletion3 = true;
 var _didIteratorError3 = false;
 var _iteratorError3 = undefined;

 try {
 var _loop = function _loop() {
 var _step3$value = _slicedToArray(_step3.value, 2),
 termIndex = _step3$value[0],
 termData = _step3$value[1];

 var _termData = _slicedToArray(termData, 7),
 ws = _termData[0],
 cs = _termData[1],
 alphas = _termData[2],
 betas = _termData[3],
 aS = _termData[4],
 bS = _termData[5],
 z = _termData[6]; // tuple


 if (verbose) {
  console.log();
  console.log('Evaluating term %i/%i : %iF%i', termIndex + 1, terms.length, aS.length, bS.length);
  console.log('powers: ws: ' + ws + ' | cs: ' + cs);
  console.log('gamma: alphas: ' + alphas + ' | betas: ' + betas);
  console.log('hyper: as: ' + aS + ' | bs: ' + bS + ' ');
  console.log('z', ctx.nstr(z));
 }

 // ===========================
 var v = hyper(aS, bS, z, opts);
 alphas.map(function (a) {
 v = _complex2.default.mul(v, math2.gamma(a));
 });
 betas.map(function (b) {
 v = _complex2.default.mul(v, math2.rgamma(b));
 });
 py.zipmap(ws, cs, function (w, c) {
 v = _complex2.default.mul(v, _complex2.default.pow(w, c));
 });
 // ===========================

 if (verbose) {
  console.log('Value: ' + v);
 }
 evaluatedTerms.push(v);
 };

 for (var _iterator3 = terms.entries()[Symbol.iterator](), _step3; !(_iteratorNormalCompletion3 = (_step3 = _iterator3.next()).done); _iteratorNormalCompletion3 = true) {
 _loop();
 }
 } catch (err) {
 _didIteratorError3 = true;
 _iteratorError3 = err;
 } finally {
 try {
 if (!_iteratorNormalCompletion3 && _iterator3.return) {
 _iterator3.return();
 }
 } finally {
 if (_didIteratorError3) {
 throw _iteratorError3;
 }
 }
 }

 if (terms.length === 1 && !perturb) {
 sumvalue = evaluatedTerms[0];
 break;
 }

 if (ctx._fixedPrecision) {
 sumvalue = _complex2.default.sum(evaluatedTerms);
 break;
 }

 sumvalue = _complex2.default.sum(evaluatedTerms);
 var termMagnitudes = evaluatedTerms.map(function (x) {
 return ctx.mag(x);
 });
 var maxMagnitude = Math.max.apply(Math, _toConsumableArray(termMagnitudes));
 var sumMagnitude = ctx.mag(sumvalue);
 var cancellation = maxMagnitude - sumMagnitude;
 if (verbose) {
  console.log('Cancellation: ' + cancellation + ' bits');
  console.log('Increased precision: ' + (ctx.prec - orig) + ' bits');
 }

 var precisionOk = cancellation < ctx.prec - orig;

 if (opts.zeroprec) {
 zeroOk = false;
 } else {
 zeroOk = maxMagnitude - ctx.prec < -zeroprec;
 }
 if (opts.infprec) {
 infOk = false;
 } else {
 infOk = maxMagnitude > infprec;
 }

 if (precisionOk && (!perturb || isNaN(cancellation))) {
 break;
 } else if (precisionOk) {
 if (perturbedReferenceValue === null) {
 hextra += 20;
 perturbedReferenceValue = sumvalue;
 continue;
 } else if (ctx.mag(sumvalue - perturbedReferenceValue) <= ctx.mag(sumvalue) - orig) {
 break;
 } else if (zeroOk) {
 sumvalue = 0;
 break;
 } else if (infOk) {
 sumvalue = Infinity;
 break;
 } else if ('hmag' in opts) {
 break;
 } else {
 hextra *= 2;
 perturbedReferenceValue = sumvalue;
 }
 } else {
 // Increase precision
 var maxDiff = Math.max(cancellation, Math.floor(orig / 2));
 var increment = Math.min(maxDiff, Math.max(extraprec, orig));
 ctx.prec += increment;
 if (verbose) {
  console.log('Must start over with increased precision');
 }
 continue;
 }
 } // while
 } finally {
 // reset precision for other funcs
 ctx.prec = orig;
 }
 return sumvalue;
}

function _checkNeedPerturb(terms, prec, discardKnownZeros) {
 var recompute = false;
 var perturb = false;
 var extraprec = 0;
 var discard = [];
 var _iteratorNormalCompletion4 = true;
 var _didIteratorError4 = false;
 var _iteratorError4 = undefined;

 try {
 for (var _iterator4 = terms.entries()[Symbol.iterator](), _step4; !(_iteratorNormalCompletion4 = (_step4 = _iterator4.next()).done); _iteratorNormalCompletion4 = true) {
 var _step4$value = _slicedToArray(_step4.value, 2),
 _termIndex = _step4$value[0],
 term = _step4$value[1];

 var _term = _slicedToArray(term, 7),
 _ws = _term[0],
 _cs = _term[1],
 _alphas = _term[2],
 _betas = _term[3],
 as = _term[4],
 bs = _term[5],
 _z2 = _term[6]; // eslint-disable-line no-unused-vars


 var haveSingularNongammaWeight = false;
 // Avoid division by zero in leading factors (TODO:
 // also check for near division by zero?)
 var _iteratorNormalCompletion5 = true;
 var _didIteratorError5 = false;
 var _iteratorError5 = undefined;

 try {
 for (var _iterator5 = _ws.entries()[Symbol.iterator](), _step5; !(_iteratorNormalCompletion5 = (_step5 = _iterator5.next()).done); _iteratorNormalCompletion5 = true) {
 var _step5$value = _slicedToArray(_step5.value, 2),
 k = _step5$value[0],
 w = _step5$value[1];

 if (_complex2.default.isZero(w)) {
 if (_complex2.default.re(_cs[k]) <= 0 && _cs[k]) {
 perturb = true;
 recompute = true;
 haveSingularNongammaWeight = true;
 }
 }
 }
 } catch (err) {
 _didIteratorError5 = true;
 _iteratorError5 = err;
 } finally {
 try {
 if (!_iteratorNormalCompletion5 && _iterator5.return) {
 _iterator5.return();
 }
 } finally {
 if (_didIteratorError5) {
 throw _iteratorError5;
 }
 }
 }

 var poleCount = [0, 0, 0];
 // Check for gamma and series poles and near-poles
 var _iteratorNormalCompletion6 = true;
 var _didIteratorError6 = false;
 var _iteratorError6 = undefined;

 try {
 for (var _iterator6 = [_alphas, _betas, bs].entries()[Symbol.iterator](), _step6; !(_iteratorNormalCompletion6 = (_step6 = _iterator6.next()).done); _iteratorNormalCompletion6 = true) {
 var _step6$value = _slicedToArray(_step6.value, 2),
 dataIndex = _step6$value[0],
 data = _step6$value[1];

 var _iteratorNormalCompletion7 = true;
 var _didIteratorError7 = false;
 var _iteratorError7 = undefined;

 try {
 for (var _iterator7 = data.entries()[Symbol.iterator](), _step7; !(_iteratorNormalCompletion7 = (_step7 = _iterator7.next()).done); _iteratorNormalCompletion7 = true) {
 var _step7$value = _slicedToArray(_step7.value, 2),
 index = _step7$value[0],
 x = _step7$value[1];

 // eslint-disable-line no-unused-vars
 var _ctx$nintDistance = ctx.nintDistance(x),
 _ctx$nintDistance2 = _slicedToArray(_ctx$nintDistance, 2),
 n = _ctx$nintDistance2[0],
 d = _ctx$nintDistance2[1];
 // Poles


 if (n > 0) {
 continue;
 }
 if (d === -Infinity) {
 // OK if we have a polynomial
 // ------------------------------
 var ok = false;
 if (dataIndex === 2) {
 var _iteratorNormalCompletion8 = true;
 var _didIteratorError8 = false;
 var _iteratorError8 = undefined;

 try {
 for (var _iterator8 = as[Symbol.iterator](), _step8; !(_iteratorNormalCompletion8 = (_step8 = _iterator8.next()).done); _iteratorNormalCompletion8 = true) {
 var u = _step8.value;

 if (ctx.isnpint(u) && u >= n) {
 ok = true;
 break;
 }
 }
 } catch (err) {
 _didIteratorError8 = true;
 _iteratorError8 = err;
 } finally {
 try {
 if (!_iteratorNormalCompletion8 && _iterator8.return) {
 _iterator8.return();
 }
 } finally {
 if (_didIteratorError8) {
 throw _iteratorError8;
 }
 }
 }
 }

 if (ok) {
 continue;
 }
 poleCount[dataIndex] += 1;
 // ------------------------------
 // perturb = recompute = True
 // return perturb, recompute, extraprec}
 } else if (d < -4) {
 extraprec += -d;
 recompute = true;
 }
 }
 } catch (err) {
 _didIteratorError7 = true;
 _iteratorError7 = err;
 } finally {
 try {
 if (!_iteratorNormalCompletion7 && _iterator7.return) {
 _iterator7.return();
 }
 } finally {
 if (_didIteratorError7) {
 throw _iteratorError7;
 }
 }
 }
 }
 } catch (err) {
 _didIteratorError6 = true;
 _iteratorError6 = err;
 } finally {
 try {
 if (!_iteratorNormalCompletion6 && _iterator6.return) {
 _iterator6.return();
 }
 } finally {
 if (_didIteratorError6) {
 throw _iteratorError6;
 }
 }
 }

 if (discardKnownZeros && poleCount[1] > poleCount[0] + poleCount[2] && !haveSingularNongammaWeight) {
 discard.append(_termIndex);
 } else if (py.sum(poleCount)) {
 perturb = true;
 recompute = true;
 }
 }
 } catch (err) {
 _didIteratorError4 = true;
 _iteratorError4 = err;
 } finally {
 try {
 if (!_iteratorNormalCompletion4 && _iterator4.return) {
 _iterator4.return();
 }
 } finally {
 if (_didIteratorError4) {
 throw _iteratorError4;
 }
 }
 }

 return [perturb, recompute, extraprec, discard];
}

function hyper(aS, bS, z, opts) {
 // """
 // Hypergeometric function, general case.
 // """
 z = _complex2.default.ensureComplex(z);
 var p = aS.length;
 var q = bS.length;
 aS = aS.map(function (a) {
 return ctx.convertParam(a);
 });
 bS = bS.map(function (b) {
 return ctx.convertParam(b);
 });

 // Reduce degree by eliminating common parameters
 if (py.get(opts, 'eliminate', true)) {
 var elimNonpositive = py.get(opts, 'eliminate_all', false);
 var i = 0;
 while (i < q && aS.length > 0) {
 var b = bS[i];
 // better way than index of to determine if b is in aS?
 var indexOfB = aS.indexOf(b);
 if (indexOfB > -1 && (elimNonpositive || !ctx.isnpint(b[0]))) {
 aS.splice(indexOfB);
 bS.splice(bS.indexOf(b));
 p -= 1;
 q -= 1;
 } else {
 i += 1;
 }
 }
 }

 // Handle special cases
 if (p === 0) {
 if (q === 1) return _hyp0f1(bS, z, opts);else if (q === 0) return _complex2.default.exp(z);
 } else if (p === 1) {
 if (q === 1) {
 var res = _hyp1f1(aS[0], bS[0], z, opts);
 return res;
 } else if (q === 2) return _hyp1f2(aS, bS, z, opts);else if (q === 0) return _hyp1f0(aS[0][0], z);
 } else if (p === 2) {
 if (q === 1) return _hyp2f1(aS, bS, z, opts);else if (q === 2) return _hyp2f2(aS, bS, z, opts);else if (q === 3) return _hyp2f3(aS, bS, z, opts);else if (q === 0) {
 var _res = _hyp2f0(aS, bS, z, opts);
 return _res;
 }
 } else if (p === q + 1) {
 return _hypq1fq(p, q, aS, bS, z, opts);
 } else if (p > q + 1 && !py.get(opts, 'force_series')) {
 return _hypBorel(p, q, aS, bS, z, opts);
 }

 var _py$unzip = py.unzip([].concat(_toConsumableArray(aS), _toConsumableArray(bS))),
 _py$unzip2 = _slicedToArray(_py$unzip, 2),
 coeffs = _py$unzip2[0],
 types = _py$unzip2[1];

 return ctx.hypsum(p, q, types, coeffs, z, opts);
}

// function hyp0f1(b, z, opts) {
// return hyper([], [b], z, opts);
// }

// function hyp1f1(a, b, z, opts) {
// return hyper([a], [b], z, opts);
// }

function hyp1f2(a1, b1, b2, z, opts) {
 return hyper([a1], [b1, b2], z, opts);
}

function hyp2f1(a, b, c, z, opts) {
 return hyper([a, b], [c], z, opts);
}

// function hyp2f2(a1, a2, b1, b2, z, opts) {
// return hyper([a1, a2], [b1, b2], z, opts);
// }

// function hyp2f3(a1, a2, b1, b2, b3, z, opts) {
// return hyper([a1, a2], [b1, b2, b3], z, opts);
// }

function hyp2f0(a, b, z, opts) {
 return hyper([a, b], [], z, opts);
}

// function hyp3f2(a1, a2, a3, b1, b2, z, opts) {
// return hyper([a1, a2, a3], [b1, b2], z, opts);
// }

// TODO
function _hyp0f1(bS, z, opts) {
 var magz = void 0,
 v = void 0;

 var _bS$ = _slicedToArray(bS[0], 2),
 b = _bS$[0],
 btype = _bS$[1];

 if (z) {
 magz = ctx.mag(z);
 } else {
 magz = 0;
 }
 if (magz >= 8 && !opts.force_series) {
 try {
 // http://functions.wolfram.com/HypergeometricFunctions/
 // Hypergeometric0F1/06/02/03/0004/
 // TODO: handle the all-real case more efficiently!
 // TODO: figure out how much precision is needed (exponential growth)
 var orig = ctx.prec;
 try {
 ctx.prec += 12 + Math.floor(magz / 2);
 var w = _complex2.default.sqrt(_complex2.default.mul(-1, z));
 var jw = _complex2.default.mul([0, 1], w);
 var minusJw = _complex2.default.mul(-1, jw);
 var u = _complex2.default.inverse(_complex2.default.mul(4, jw));
 var minusU = _complex2.default.mul(-1, u);
 var c = _complex2.default.sub(0.5, b);
 var bMinusHalf = _complex2.default.sub(b, 0.5);
 var oneAndHalfMinusB = _complex2.default.sub(1.5, b);
 var capitalE = _complex2.default.exp(_complex2.default.mul(2, jw));
 var hyp2f0Res = hyp2f0(bMinusHalf, oneAndHalfMinusB, minusU, { 'force_series': true });
 var prefactor = _complex2.default.div(_complex2.default.pow(minusJw, c), capitalE);
 var capitalH1 = _complex2.default.mul(prefactor, hyp2f0Res);
 hyp2f0Res = hyp2f0(bMinusHalf, oneAndHalfMinusB, u, { 'force_series': true });
 prefactor = _complex2.default.mul(_complex2.default.pow(jw, c), capitalE);
 var capitalH2 = _complex2.default.mul(prefactor, hyp2f0Res);
 var hSum = _complex2.default.add(capitalH1, capitalH2);
 var quotient = _complex2.default.div(math2.gamma(b), 2 * Math.sqrt(Math.PI));
 v = _complex2.default.mul(quotient, hSum);
 } finally {
 ctx.prec = orig;
 }
 // if b and z are real
 if (ctx._isRealType(b) && ctx._isRealType(z)) {
 v = _complex2.default.re(v);
 }
 // below was return +v - in python this means carry context settings like rounding out with the decimal.
 return v;
 } catch (e) {
 if (e.name !== 'NoConvergence') {
 throw e;
 }
 }
 }
 return ctx.hypsum(0, 1, [btype], [b], z, opts);
};

function _hyp1f1(aS, bS, z, opts) {
 var _aS = _slicedToArray(aS, 2),
 a = _aS[0],
 atype = _aS[1];

 var _bS = _slicedToArray(bS, 2),
 b = _bS[0],
 btype = _bS[1];

 if (_complex2.default.isZero(z)) {
 return _complex2.default.add(1, z);
 }
 var magz = ctx.mag(z);
 var v = void 0;
 if (magz >= 7 && !(ctx.isint(a) && _complex2.default.re(a) <= 0)) {
 if (!isFinite(z)) {
 if (mpFunctions.sign(a) === mpFunctions.sign(b) === mpFunctions.sign(z) === 1) {
 return Infinity;
 }
 return _complex2.default.mul(NaN, z);
 }
 try {
 try {
 ctx.prec += magz;
 var sector = _complex2.default.im(z) < 0;
 var h = function h(a, b) {
 var negA = _complex2.default.mul(-1, a);
 var aMinusB = _complex2.default.sub(a, b);
 var bMinusA = _complex2.default.sub(b, a);
 var capitalE = void 0;
 if (sector) {
 capitalE = ctx.expjpi(negA);
 } else {
 capitalE = ctx.expjpi(a);
 }
 var rz = _complex2.default.inverse(z);
 var capitalT1 = [[capitalE, z], [1, negA], [b], [bMinusA], [a, _complex2.default.add(1, aMinusB)], [], _complex2.default.mul(-1, rz)];
 var capitalT2 = [[_complex2.default.exp(z), z], [1, aMinusB], [b], [a], [bMinusA, 1 - a], [], rz];
 return [capitalT1, capitalT2];
 };
 v = hypercomb(h, [a, b], { 'force_series': true });
 if (ctx._is_real_type(a) && ctx._is_real_type(b) && ctx._is_real_type(z)) {
 v = _complex2.default.re(v);
 }
 return v;
 } catch (e) {
 if (e.name !== 'NoConvergence') {
 throw e;
 }
 }
 } finally {
 ctx.prec -= magz;
 }
 }

 v = ctx.hypsum(1, 1, (atype, btype), [a, b], z, opts);
 return v;
};

function _hyp2f1Gosper(a, b, c, z) {
 var opts = arguments.length > 4 && arguments[4] !== undefined ? arguments[4] : {};

 var f1 = void 0;
 // Use Gosper's recurrence
 // See http://www.math.utexas.edu/pipermail/maxima/2006/000126.html
 var _a = a,
 _b = b,
 _c = c,
 _z = z;

 var orig = ctx.prec;
 var maxprec = py.get(opts, 'maxprec', 100 * orig);
 var extra = 10;
 while (1) {
 ctx.prec = orig + extra;
 a = ctx.convert(_a);
 b = ctx.convert(_b);
 c = ctx.convert(_c);
 z = ctx.convert(_z);
 var d = 0;
 var e = 1;
 var f = 0;
 var k = 0;
 // Common subexpression elimination, unfortunately making
 // things a bit unreadable. The formula is quite messy to begin
 // with, though...
 var abz = _complex2.default.mul(_complex2.default.mul(a, b), z);
 var ch = _complex2.default.mul(c, 0.5);
 var c1h = _complex2.default.mul(_complex2.default.add(c, 1), 0.5);
 var nz = _complex2.default.sub(1, z);
 var g = _complex2.default.div(z, nz);
 var abg = _complex2.default.mul(_complex2.default.mul(a, b), g);
 var cba = _complex2.default.sub(_complex2.default.sub(c, b), a);
 var z2 = _complex2.default.sub(z, 2);
 var tol = -ctx.prec - 10;
 var maxmag = -Infinity;
 while (1) {
 var kch = _complex2.default.add(k, ch);
 // kakbz = (k+a)*(k+b)*z / (4*(k+1)*kch*(k+c1h));
 var numerator = _complex2.default.prod([_complex2.default.add(k, a), _complex2.default.add(k, b), z]);
 var denominator = _complex2.default.prod([4, _complex2.default.add(k, 1), kch, _complex2.default.add(k, c1h)]);
 var kakbz = _complex2.default.div(numerator, denominator);
 // d1 = kakbz*(e-(k+cba)*d*g)
 var kcbadg = _complex2.default.prod([_complex2.default.add(k, cba), d, g]);
 var d1 = _complex2.default.mul(kakbz, _complex2.default.sub(e, kcbadg));
 // e1 = kakbz*(d*abg+(k+c)*e);
 var dabgkce = _complex2.default.add(_complex2.default.mul(d, abg), _complex2.default.mul(_complex2.default.add(k, c), e));
 var e1 = _complex2.default.mul(kakbz, dabgkce);
 // ft = d*(k*(cba*z+k*z2-c)-abz)/(2*kch*nz);
 var cbazkz2 = _complex2.default.add(_complex2.default.mul(cba, z), _complex2.default.mul(k, z2));
 var cbazkz2c = _complex2.default.sub(cbazkz2, c);
 var kcbazkz2cabz = _complex2.default.sub(_complex2.default.mul(k, cbazkz2c), abz);
 numerator = _complex2.default.mul(d, kcbazkz2cabz);
 denominator = _complex2.default.prod([2, kch, nz]);
 var ft = _complex2.default.div(numerator, denominator);
 f1 = _complex2.default.sub(_complex2.default.add(f, e), ft);
 maxmag = Math.max(maxmag, ctx.mag(f1));
 if (ctx.mag(_complex2.default.sub(f1, f)) < tol) {
 break;
 }
 d = d1;
 e = e1;
 f = f1;

 k = _complex2.default.add(k, 1);
 }
 var cancellation = maxmag - ctx.mag(f1);
 if (cancellation < extra) {
 break;
 } else {
 extra += cancellation;
 if (extra > maxprec) {
 throw new ctx.NoConvergence();
 }
 }
 }
 return f1;
}

function _hyp1f2(aS, bS, z, opts) {
 var _aS$ = _slicedToArray(aS[0], 2),
 a1 = _aS$[0],
 a1type = _aS$[1];

 var _bS2 = _slicedToArray(bS, 2),
 _bS2$ = _slicedToArray(_bS2[0], 2),
 b1 = _bS2$[0],
 b1type = _bS2$[1],
 _bS2$2 = _slicedToArray(_bS2[1], 2),
 b2 = _bS2$2[0],
 b2type = _bS2$2[1];

 var absz = _complex2.default.abs(z);
 var magz = ctx.mag(z);
 var orig = ctx.prec;

 // Asymptotic expansion is ~ exp(sqrt(z))
 var asympExtraprec = z && Math.floor(magz / 2);

 // # Asymptotic series is in terms of 3F0
 // set to true for the moment, normally false for some reason.
 var canUseAsymptotic = !py.get(opts, 'force_series') && ctx.mag(absz) > 19 && Math.sqrt(absz) > 1.5 * orig; // # and \
 // # ctx._hyp_check_convergence([a1, a1-b1+1, a1-b2+1], [],
 // # 1/absz, orig+40+asympExtraprec)

 // # TODO: much of the following could be shared with 2F3 instead of
 // # copypasted
 if (canUseAsymptotic) {
 // #print "using asymp"
 try {
 try {
 ctx.prec += asympExtraprec;
 // # http://functions.wolfram.com/HypergeometricFunctions/
 // # Hypergeometric1F2/06/02/03/
 var h = function h(a1, b1, b2) {
 var a1Sqrd = _complex2.default.pow(a1, 2);
 var b1Plusb2 = _complex2.default.add(b1, b2);
 var a1Minusb1Minusb2 = _complex2.default.sub(a1, b1Plusb2);
 var capitalX = _complex2.default.mul(0.5, _complex2.default.add(a1Minusb1Minusb2, 0.5));
 var c = [];
 c[0] = 1;
 // c1
 var b1Plusb2Minus2 = _complex2.default.sub(b1Plusb2, 2);
 var threeA1Plus = _complex2.default.add(_complex2.default.mul(3, a1), b1Plusb2Minus2);
 var b1b2 = _complex2.default.mul(b1, b2);
 var prod = _complex2.default.prod([0.25, threeA1Plus, a1Minusb1Minusb2]);
 var bigTerm = _complex2.default.sum(prod, b1b2, -3 / 16);
 c[1] = _complex2.default.mul(2, bigTerm);

 // c2
 var nEightA1Sqr = _complex2.default.mul(-8, a1Sqrd);
 var sum = _complex2.default.sum([nEightA1Sqr, _complex2.default.mul(11, a1), b1Plusb2Minus2]);
 var twoA1minus3 = _complex2.default.sub(_complex2.default.mul(2, a1), 3);
 var capitalA = _complex2.default.prod([-16, twoA1minus3, b1b2]);
 var capitalB = _complex2.default.prod([4, a1Minusb1Minusb2, sum]);
 var thing3 = _complex2.default.sum(capitalA, capitalB, -3);
 var twoBigTermSqrd = _complex2.default.mul(2, _complex2.default.pow(bigTerm, 2));
 c[2] = _complex2.default.add(twoBigTermSqrd, _complex2.default.mul(1 / 16, thing3));
 var s1 = 0;
 var s2 = 0;
 var k = 0;
 var tprev = 0;

 while (1) {
 if (c[k] === undefined) {
 var _sum = _complex2.default.sum([_complex2.default.mul(-6, a1), _complex2.default.mul(2, b1), _complex2.default.mul(2, b2), -4]);
 var kSum = _complex2.default.mul(k, _sum);
 var threeA1Sqrd = _complex2.default.mul(3, a1Sqrd);
 var threekSqrd = _complex2.default.mul(3, _complex2.default.pow(k, 2));
 var uu1 = _complex2.default.sum([threekSqrd, kSum, threeA1Sqrd, _complex2.default.mul(-1, _complex2.default.pow(b1Plusb2, 2)), _complex2.default.prod([-1, 2, a1, b1Plusb2Minus2]), 0.25]);
 var minusA1 = _complex2.default.mul(-1, a1);
 var minusB1 = _complex2.default.mul(-1, b1);
 var minusB2 = _complex2.default.mul(-1, b2);
 var uu2 = _complex2.default.prod([_complex2.default.sum([k, minusA1, b1, minusB2, -0.5]), _complex2.default.sum([k, minusA1, minusB1, b2, -0.5]), _complex2.default.sum([k, minusA1, b1, b2, -5 / 2])]);
 var oneOver2k = _complex2.default.inverse(_complex2.default.mul(2, k));
 var diff = _complex2.default.sub(_complex2.default.mul(uu1, c[k - 1]), _complex2.default.mul(uu2, c[k - 2]));
 c[k] = _complex2.default.mul(oneOver2k, diff);
 }
 var minusHalfK = _complex2.default.mul(-0.5, k);
 var twoToTheMinusK = _complex2.default.pow(2, _complex2.default.mul(-1, k));
 var minusZ = _complex2.default.mul(-1, z);
 var w = _complex2.default.mul(c[k], _complex2.default.pow(minusZ, minusHalfK));
 var t1 = _complex2.default.prod([_complex2.default.pow([0, -1], k), twoToTheMinusK, w]);
 var t2 = _complex2.default.prod([_complex2.default.pow(ctx.j, k), twoToTheMinusK, w]);
 if (_complex2.default.abs(t1) < 0.1 * ctx.eps) {
 // #print "Convergence :)"
 break;
 }
 // # Quit if the series doesn't converge quickly enough
 if (k > 5 && _complex2.default.abs(tprev) / _complex2.default.abs(t1) < 1.5) {
 // #print "No convergence :("
 throw new ctx.NoConvergence();
 }
 s1 = _complex2.default.add(s1, t1);
 s2 = _complex2.default.add(s2, t2);
 tprev = t1;
 k += 1;
 }
 var capitalS = ctx.expj(Math.PI * capitalX + 2 * _complex2.default.sqrt(-z)) * s1 + ctx.expj(-(Math.PI * capitalX + 2 * _complex2.default.sqrt(-z))) * s2;
 var capitalT1 = [[0.5 * capitalS, Math.PI, -z], [1, -0.5, capitalX], [b1, b2], [a1], [], [], 0];
 var capitalT2 = [[-z], [-a1], [b1, b2], [b1 - a1, b2 - a1], [a1, a1 - b1 + 1, a1 - b2 + 1], [], 1 / z];
 return [capitalT1, capitalT2];
 };

 var _v = hypercomb(h, [a1, b1, b2], { force_series: true, maxterms: 4 * ctx.prec });
 if ([a1, b1, b2, z].every(ctx._isRealType)) {
 _v = _complex2.default.re(_v);
 }
 return _v;
 } catch (e) {
 if (e.name === 'NoConvergence') {
 // pass // python statement for do nothing.
 } else {
 // only catch NoConvergence errors
 throw e;
 }
 }
 } finally {
 ctx.prec = orig;
 }
 }
 // #print "not using asymp"
 var res = ctx.hypsum(1, 2, (a1type, b1type, b2type), [a1, b1, b2], z, opts);
 return res;
};
function _hyp1f0() {
 throw new Error('Not Implemented');
};
function _hyp2f1(aS, bS, z) {
 var opts = arguments.length > 3 && arguments[3] !== undefined ? arguments[3] : {};

 var h = void 0,
 v = void 0,
 t = void 0,
 ca = void 0,
 cb = void 0,
 rz = void 0,
 capitalT1 = void 0,
 capitalT2 = void 0;

 var _aS2 = _slicedToArray(aS, 2),
 _aS2$ = _slicedToArray(_aS2[0], 2),
 a = _aS2$[0],
 atype = _aS2$[1],
 _aS2$2 = _slicedToArray(_aS2[1], 2),
 b = _aS2$2[0],
 btype = _aS2$2[1];

 var _bS3 = _slicedToArray(bS, 2),
 c = _bS3[0],
 ctype = _bS3[1];

 if (_complex2.default.equals(z, 1)) {
 // TODO: the following logic can be simplified
 var convergent = _complex2.default.re(_complex2.default.sub(c, _complex2.default.add(a, b))) > 0;
 var finite = ctx.isint(a) && a <= 0 || ctx.isint(b) && b <= 0;
 var condition2 = !(ctx.isint(a) && c <= a <= 0 || ctx.isint(b) && c <= b <= 0);
 var zerodiv = ctx.isint(c) && c <= 0 && condition2;
 // print "bz", a, b, c, z, convergent, finite, zerodiv
 // Gauss's theorem gives the value if convergent
 if ((convergent || finite) && !zerodiv) {
 var _ca = _complex2.default.sub(c, a);
 return mpFactorials.gammaprod([c, _complex2.default.minus(_ca, b)], [_ca, _complex2.default.sub(c, b)], true);
 }
 // Otherwise, there is a pole and we take the
 // sign to be that when approaching from below
 // XXX: this evaluation is not necessarily correct in all cases
 return _complex2.default.mul(hyp2f1(a, b, c, 1 - ctx.eps * 2), Infinity);
 }
 // Equal to 1 (first term), unless there is a subsequent
 // division by zero
 if (!_complex2.default.isZero(z)) {
 // Division by zero but power of z is higher than
 // first order so cancels
 if (!_complex2.default.isZero(c) || _complex2.default.isZero(a) || _complex2.default.isZero(b)) {
 return _complex2.default.add(1, z);
 }
 // Indeterminate
 return NaN;
 }

 // Hit zero denominator unless numerator goes to 0 first
 if (ctx.isint(c) && c <= 0) {
 if (ctx.isint(a) && c <= a <= 0 || ctx.isint(b) && c <= b <= 0) {
 // pass; // Python statement to move on.
 } else {
 // Pole in series
 return Infinity;
 }
 }
 var absz = _complex2.default.abs(z);

 // Fast case: standard series converges rapidly,
 // possibly in finitely many terms
 if (absz <= 0.8 || ctx.isint(a) && a <= 0 && a >= -1000 || ctx.isint(b) && b <= 0 && b >= -1000) {
 return ctx.hypsum(2, 1, [atype, btype, ctype], [a, b, c], z, opts);
 }
 var orig = ctx.prec;
 try {
 ctx.prec += 10;
 // Use 1/z transformation
 if (absz >= 1.3) {
 h = function h(a, b) {
 var t = _complex2.default.sub(1, c);
 var ab = _complex2.default.sub(a, b);
 var rz = _complex2.default.div(1, z);
 var capitalT1 = ([_complex2.default.mul(-1, z)], [_complex2.default.mul(-1, a)], [c, _complex2.default.mul(-1, ab)], [b, _complex2.default.sub(c, a)], [a, _complex2.default.add(t, a)], [_complex2.default.add(1, ab)], rz);
 var capitalT2 = ([_complex2.default.mul(-1, z)], [_complex2.default.mul(-1, b)], [c, ab], [a, _complex2.default.sub(c, b)], [b, _complex2.default.add(t, b)], [_complex2.default.sub(1, ab)], rz);
 return [capitalT1, capitalT2];
 };
 v = hypercomb(h, [a, b], opts);

 // Use 1-z transformation
 } else if (_complex2.default.abs(_complex2.default.sub(1, z)) <= 0.75) {
 h = function h(a, b) {
 t = _complex2.default.sub(c, _complex2.default.add(a, b));
 ca = _complex2.default.sub(c, a);
 cb = _complex2.default.sub(c, b);
 rz = _complex2.default.sub(1, z);
 capitalT1 = [[], [], [c, t], [ca, cb], [a, b], [_complex2.default.sub(1, t)], rz];
 capitalT2 = [[rz], [t], [c, _complex2.default.mul(-1, t)], [a, b], [ca, cb], [_complex2.default.add(1, t)], rz];
 return [capitalT1, capitalT2];
 };
 v = hypercomb(h, [a, b], opts);

 // Use z/(z-1) transformation
 } else if (_complex2.default.abs(_complex2.default.div(z, _complex2.default.sub(z, 1))) <= 0.75) {
 var h2f1 = hyp2f1(a, _complex2.default.sub(c, b), c, _complex2.default.div(z, _complex2.default.sub(z, 1)));
 v = _complex2.default.div(h2f1, _complex2.default.pow(_complex2.default.sub(1, z), a));
 } else {
 // Remaining part of unit circle
 v = _hyp2f1Gosper(a, b, c, z, opts);
 }
 } finally {
 ctx.prec = orig;
 }
 return v;
};

function _hyp2f2() {
 throw new Error('Not Implemented');
};
function _hyp2f3() {
 throw new Error('Not Implemented');
};
function _hyp2f0(aS, bS, z) {
 var opts = arguments.length > 3 && arguments[3] !== undefined ? arguments[3] : {};

 var _aS3 = _slicedToArray(aS, 2),
 _aS3$ = _slicedToArray(_aS3[0], 2),
 a = _aS3$[0],
 atype = _aS3$[1],
 _aS3$2 = _slicedToArray(_aS3[1], 2),
 b = _aS3$2[0],
 btype = _aS3$2[1];
 // We want to try aggressively to use the asymptotic expansion,
 // and fall back only when absolutely necessary


 try {
 var optsb = JSON.parse(JSON.stringify(opts)); // TODO: do we really need a copy? copying is sorta rough in JS
 optsb['maxterms'] = py.get(optsb, 'maxterms', ctx.prec);
 return ctx.hypsum(2, 0, [atype, btype], [a, b], z, optsb);
 } catch (e) {
 if (e.name === 'NoConvergence') {
 if (py.get(opts, 'force_series')) {
 throw e;
 }
 } else {
 // only catch NoConvergence errors
 throw e;
 }
 }

 // THOUGHTS: you have the perturb of a so handle the subtraction better down below?
 var h = function h(a, b) {
 var opts = arguments.length > 2 && arguments[2] !== undefined ? arguments[2] : {};

 // HACK = manually handle perturbations
 // sin(pi*(b + h)) = sin(pi*b)cos(pi*h) + cos(pi * b)sin(pi*h)
 //
 // b = 0, 2, 4,
 // sin(pi*h)
 //
 // b = 1,3,5
 // -sin(pi*h)
 //
 // b = 1/2, 5/2, 9/2
 // cos(pi*h)
 //
 // b = 3/2, 7/2, 11/2
 // -cos(pi*h)
 var w = void 0;
 if (opts.perturbParams) {
 switch ((b * 2 + 1) % 4) {
 case 1:
 // 0, 2, 4
 w = Math.sin(Math.PI * opts.perturbParams[1]);
 break;
 case 2:
 // 1/2, 5/2, 9/2
 w = Math.cos(Math.PI * opts.perturbParams[1]);
 break;
 case 3:
 // 1,3,5
 w = -Math.sin(Math.PI * opts.perturbParams[1]);
 break;
 case 0:
 // 3/2, 7/2, 11/2
 w = -Math.cos(Math.PI * opts.perturbParams[1]);
 break;
 default:
  // console.log('sinpi was perturbed but b was not half integer');
  w = math2.sinpi(b);
 }
 } else {
 w = math2.sinpi(b);
 }
 var rz = _complex2.default.mul(-1, _complex2.default.inverse(z));
 var aMinusB = _complex2.default.sub(a, b);
 var twoMinusB = _complex2.default.sub(2, b);
 // ws, cs, alphas, betas, as, bs, z
 var capitalT1 = [[Math.PI, w, rz], [1, -1, a], [], [_complex2.default.add(aMinusB, 1), b], [a], [b], rz];
 var capitalT2 = [[-Math.PI, w, rz], [1, -1, _complex2.default.add(1, aMinusB)], [], [a, twoMinusB], [_complex2.default.add(aMinusB, 1)], [twoMinusB], rz];
 return [capitalT1, capitalT2];
 };

 var aMinusB = _complex2.default.sub(a, b);
 return hypercomb(h, [a, _complex2.default.add(1, aMinusB)], opts);
};

function _hypq1fq() {
 throw new Error('Not Implemented');
};

// TODO - this was abandoned near completion because we realized we don't need it for
// the use case we're starting with. It's only lacking the quad function.
function _hypBorel(p, q, aS, bS, z, opts) {
 throw new Error('Not Fully Implemented / Tested');
 // if (aS) {
 // [aS, aTypes] = pyp.unzip(aS)
 // } else {
 // aS = []
 // aTypes = []
 // }
 // if (bS) {
 // [bS, bTypes] = py.unzip(bS)
 // } else {
 // bS = []
 // bTypes = []
 // }
 //
 // opts['maxterms'] = py.get(opts, 'maxterms', ctx.prec)
 // try {
 // return ctx.hypsum(p, q, aTypes + bTypes, aS + bS, z, opts)
 // } catch (e) {
 // if (e.name === 'NoConvergence') {
 // // pass // python statement for do nothing.
 // } else {
 // // only catch NoConvergence errors
 // throw e
 // }
 // }
 // prec = ctx.prec
 // try {
 // tol = opts.asymp_tol || (ctx.eps / 4)
 // ctx.prec += 10
 // function term (k, cache = { 0: 1 }) {
 // if (k in cache) {
 // return cache[k]
 // }
 // t = term(k - 1)
 // aS.map((a) => t *= a + (k - 1))
 // bS.map((b) => t /= b + (k - 1))
 // t *= z
 // t /= k
 // cache[k] = t
 // return t
 // }
 // s = 1
 // for (k of py.range(1, ctx.prec)) {
 // t = term(k)
 // s += t
 // if (Math.abs(t) <= tol) {
 // return s
 // }
 // }
 // } finally {
 // ctx.prec = prec
 // }
 // if (p <= q + 3) {
 // if (!opts.contour) {
 // if (Complex.arg(z) < 0.25) {
 // u = Complex.div(z, max(1, Complex.abs(z)))
 // if (Complex.arg(z) >= 0) {
 // contour = [0, [0, 2], Complex.div([2, 2], u), Complex.div(2, u), Infinity]
 // } else {
 // contour = [0, [0, -2], Complex.div([2, -2], u), Complex.div(2, u), Infinity]
 // }
 // // contour = [0, 2j/z, 2/z, ctx.inf]
 // // contour = [0, 2j, 2/z, ctx.inf]
 // // contour = [0, 2j, ctx.inf]}
 // } else {
 // contour = [0, Infinity]
 // }
 // }
 // quadOpts = opts.quadOpts || {}
 // function g (t) {
 // eToTheMinusT = Complex.exp(-t)
 // hyperVal = hyper(aS, [...bS, 1], Complex.mul(t, z))
 // return Complex.mul(eToTheMinusT, hyperVal)
 // }
 // quadOpts['error'] = true;
 // [capitalI, err] = mpmath.quad(g, contour, quadOpts) // todo-missingfunc quad
 // if (err <= abs(capitalI) * ctx.eps * 8) {
 // return capitalI
 // }
 // }
 // throw new ctx.UserException('_hyp_borel failed to converge')
};

exports.hypercomb = hypercomb;
exports.hyper = hyper;
exports.hyp2f0 = hyp2f0;
exports.hyp1f2 = hyp1f2;
},{"../../../utils/complex.js":89,"../../../utils/pythonHelpers.js":94,"../../math2/math2.js":83,"../ctx.js":84,"./factorials.js":86,"./functions.js":87}],89:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});

var _slicedToArray = function () { function sliceIterator(arr, i) { var _arr = []; var _n = true; var _d = false; var _e = undefined; try { for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) { _arr.push(_s.value); if (i && _arr.length === i) break; } } catch (err) { _d = true; _e = err; } finally { try { if (!_n && _i["return"]) _i["return"](); } finally { if (_d) throw _e; } } return _arr; } return function (arr, i) { if (Array.isArray(arr)) { return arr; } else if (Symbol.iterator in Object(arr)) { return sliceIterator(arr, i); } else { throw new TypeError("Invalid attempt to destructure non-iterable instance"); } }; }();

var _createClass = function () { function defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } } return function (Constructor, protoProps, staticProps) { if (protoProps) defineProperties(Constructor.prototype, protoProps); if (staticProps) defineProperties(Constructor, staticProps); return Constructor; }; }();

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

// List of sources:
//
// ComplexJS: https://github.com/infusion/Complex.js

// Thoughts:
// Inline functions instead of calling class?
var Complex = function () {
 function Complex() {
 _classCallCheck(this, Complex);
 }

 _createClass(Complex, null, [{
 key: "equals",
 value: function equals(a, b) {
 var precision = arguments.length > 2 && arguments[2] !== undefined ? arguments[2] : 1e-14;

 a = Complex.ensureComplex(a);
 b = Complex.ensureComplex(b);
 return Math.abs(a[0] - b[0]) <= precision && Math.abs(a[1] - b[1]) <= precision;
 }
 }, {
 key: "abs",
 value: function abs(a) {
 if (Complex.isReal(a)) {
 return Math.abs(Complex.re(a));
 }
 return Math.sqrt(a[0] ** 2 + a[1] ** 2);
 }

 // angle

 }, {
 key: "arg",
 value: function arg(a) {
 a = Complex.ensureComplex(a);
 return Math.atan2(a[1], a[0]);
 }
 }, {
 key: "add",
 value: function add(a, b) {
 if (Complex.isReal(a) && Complex.isReal(b)) {
 return Complex.re(a) + Complex.re(b);
 }
 a = Complex.ensureComplex(a);
 b = Complex.ensureComplex(b);
 return [a[0] + b[0], a[1] + b[1]];
 }
 }, {
 key: "sub",
 value: function sub(a, b) {
 if (Complex.isReal(a) && Complex.isReal(b)) {
 return Complex.re(a) - Complex.re(b);
 }
 a = Complex.ensureComplex(a);
 b = Complex.ensureComplex(b);
 return [a[0] - b[0], a[1] - b[1]];
 }
 }, {
 key: "mul",
 value: function mul(a, b) {
 if (Complex.isReal(a) && Complex.isReal(b)) {
 return Complex.re(a) * Complex.re(b);
 }
 a = Complex.ensureComplex(a);
 b = Complex.ensureComplex(b);
 var realPart = a[0] * b[0] - a[1] * b[1];
 var imagPart = a[0] * b[1] + a[1] * b[0];
 return [realPart, imagPart];
 }
 }, {
 key: "div",
 value: function div(a, b) {
 if (Complex.isReal(a) && Complex.isReal(b)) {
 return Complex.re(a) / Complex.re(b);
 }
 a = Complex.ensureComplex(a);
 b = Complex.ensureComplex(b);
 var denominator = b[0] ** 2 + b[1] ** 2;
 var realPart = (a[0] * b[0] + a[1] * b[1]) / denominator;
 var imagPart = (a[1] * b[0] - a[0] * b[1]) / denominator;
 return [realPart, imagPart];
 }
 }, {
 key: "sqrt",
 value: function sqrt(a) {
 if (Complex.isReal(a) && a >= 0) {
 return Math.sqrt(Complex.re(a));
 }
 if (Complex.isReal(a) && a < 0) {
 return [0, Math.sqrt(-Complex.re(a))];
 }
 a = Complex.ensureComplex(a);
 var mod = Math.sqrt(a[0] ** 2 + a[1] ** 2);
 var realPart = Math.sqrt((mod + a[0]) / 2);
 var sign = Math.sign(a[1]);
 // If imaginary part is 0 then use positive
 if (sign === 0) {
 sign = 1;
 }
 var imagPart = sign * Math.sqrt((mod - a[0]) / 2);
 return [realPart, imagPart];
 }
 }, {
 key: "exp",
 value: function exp(a) {
 if (Complex.isReal(a)) {
 return Math.exp(Complex.re(a));
 }
 var expA = Math.exp(a[0]);
 var realPart = expA * Math.cos(a[1]);
 var imagPart = expA * Math.sin(a[1]);
 return [realPart, imagPart];
 }

 /* Notes from complex.js for Complex.pow
 *
 * a + bi ^ c + di = (a + bi)^(c + di)
 * = exp((c + di) * log(a + bi)
 * = pow(a^2 + b^2, (c + di) / 2) * exp(i(c + di)atan2(b, a))
 * =>...
 * Re = (pow(a^2 + b^2, c / 2) * exp(-d * atan2(b, a))) * cos(d * log(a^2 + b^2) / 2 + c * atan2(b, a))
 * Im = (pow(a^2 + b^2, c / 2) * exp(-d * atan2(b, a))) * sin(d * log(a^2 + b^2) / 2 + c * atan2(b, a))
 *
 * =>...
 * Re = exp(c * log(sqrt(a^2 + b^2)) - d * atan2(b, a)) * cos(d * log(sqrt(a^2 + b^2)) + c * atan2(b, a))
 * Im = exp(c * log(sqrt(a^2 + b^2)) - d * atan2(b, a)) * sin(d * log(sqrt(a^2 + b^2)) + c * atan2(b, a))
 *
 * =>
 * Re = exp(c * logsq2 - d * arg(z_1)) * cos(d * logsq2 + c * arg(z_1))
 * Im = exp(c * logsq2 - d * arg(z_1)) * sin(d * logsq2 + c * arg(z_1))
 *
 */

 }, {
 key: "pow",
 value: function pow(a, b) {
 if (Complex.isReal(a) && Complex.isReal(b)) {
 var result = Math.pow(Complex.re(a), Complex.re(b));
 if (!isNaN(result)) {
 return result;
 }
 }
 a = Complex.ensureComplex(a);
 b = Complex.ensureComplex(b);
 var atan2 = Math.atan2(a[1], a[0]);
 var logsqr = Math.log(Math.sqrt(a[0] ** 2 + a[1] ** 2));

 var preR = Math.exp(b[0] * logsqr - b[1] * atan2);
 var preRi = b[1] * logsqr + b[0] * atan2;

 var realPart = preR * Math.cos(preRi);
 var imagPart = preR * Math.sin(preRi);

 return [realPart, imagPart];
 }
 }, {
 key: "log",
 value: function log(a) {
 if (Complex.isReal(a)) {
 return Math.log(Complex.re(a));
 }
 var realPart = Math.log(Math.sqrt(a[0] ** 2 + a[1] ** 2));
 var imagPart = Math.atan2(a[1], a[0]);
 return [realPart, imagPart];
 }
 }, {
 key: "sin",
 value: function sin(a) {
 if (Complex.isReal(a)) {
 return Math.sin(Complex.re(a));
 }
 var realPart = Math.sin(a[0]) * Math.cosh(a[1]);
 var imagPart = Math.cos(a[0]) * Math.sinh(a[1]);
 return [realPart, imagPart];
 }
 }, {
 key: "cos",
 value: function cos(a) {
 if (Complex.isReal(a)) {
 return Math.cos(Complex.re(a));
 }
 var realPart = Math.cos(a[0]) * Math.cosh(a[1]);
 var imagPart = -Math.sin(a[0]) * Math.sinh(a[1]);
 return [realPart, imagPart];
 }
 }, {
 key: "tan",
 value: function tan(a) {
 if (Complex.isReal(a)) {
 return Math.tan(Complex.re(a));
 }
 var denominator = Math.cos(a[0] * 2) + Math.cosh(a[1] * 2);
 var realPart = Math.sin(a[0] * 2) / denominator;
 var imagPart = Math.sinh(a[1] * 2) / denominator;
 return [realPart, imagPart];
 }
 }, {
 key: "cot",
 value: function cot(a) {
 a = Complex.ensureComplex(a);
 var denominator = Math.cos(a[0] * 2) - Math.cosh(a[1] * 2);
 var realPart = -Math.sin(a[0] * 2) / denominator;
 var imagPart = Math.sinh(a[1] * 2) / denominator;
 return [realPart, imagPart];
 }
 }, {
 key: "sec",
 value: function sec(a) {
 a = Complex.ensureComplex(a);
 var denominator = 0.5 * Math.cosh(a[1] * 2) + 0.5 * Math.cos(a[0] * 2);
 var realPart = Math.cos(a[0]) * Math.cosh(a[1]) / denominator;
 var imagPart = Math.sin(a[0]) * Math.sinh(a[1]) / denominator;
 return [realPart, imagPart];
 }
 }, {
 key: "csc",
 value: function csc(a) {
 a = Complex.ensureComplex(a);
 var denominator = 0.5 * Math.cosh(a[1] * 2) - 0.5 * Math.cos(a[0] * 2);
 var realPart = Math.sin(a[0]) * Math.cosh(a[1]) / denominator;
 var imagPart = -Math.cos(a[0]) * Math.sinh(a[1]) / denominator;
 return [realPart, imagPart];
 }
 }, {
 key: "asin",
 value: function asin(a) {
 if (Complex.isReal(a)) {
 return Math.asin(Complex.re(a));
 }
 var realPart = a[1] ** 2 - a[0] ** 2 + 1;
 var imagPart = -2 * a[0] * a[1];

 var _Complex$sqrt = Complex.sqrt([realPart, imagPart]);

 var _Complex$sqrt2 = _slicedToArray(_Complex$sqrt, 2);

 realPart = _Complex$sqrt2[0];
 imagPart = _Complex$sqrt2[1];


 realPart = realPart - a[1];
 imagPart = imagPart + a[0];

 var _Complex$log = Complex.log([realPart, imagPart]);

 var _Complex$log2 = _slicedToArray(_Complex$log, 2);

 realPart = _Complex$log2[0];
 imagPart = _Complex$log2[1];


 return [imagPart, -realPart];
 }
 }, {
 key: "acos",
 value: function acos(a) {
 if (Complex.isReal(a)) {
 return Math.acos(Complex.re(a));
 }
 var realPart = void 0,
 imagPart = void 0;

 var _Complex$asin = Complex.asin(a);

 var _Complex$asin2 = _slicedToArray(_Complex$asin, 2);

 realPart = _Complex$asin2[0];
 imagPart = _Complex$asin2[1];

 return [Math.PI / 2 - realPart, -imagPart];
 }
 }, {
 key: "atan",
 value: function atan(a) {
 if (Complex.isReal(a)) {
 return Math.atan(Complex.re(a));
 }
 var denominator = a[0] ** 2 + (a[1] - 1) ** 2;

 var realPart = (1 - a[1] ** 2 - a[0] ** 2) / denominator;
 var imagPart = -2 * a[0] / denominator;

 var _Complex$log3 = Complex.log([realPart, imagPart]);

 var _Complex$log4 = _slicedToArray(_Complex$log3, 2);

 realPart = _Complex$log4[0];
 imagPart = _Complex$log4[1];


 return [imagPart * -0.5, realPart * 0.5];
 }
 }, {
 key: "acot",
 value: function acot(a) {
 a = Complex.ensureComplex(a);
 var denominator = a[0] ** 2 + a[1] ** 2;

 var realPart = a[0] / denominator;
 var imagPart = -a[1] / denominator;

 return Complex.atan([realPart, imagPart]);
 }
 }, {
 key: "asec",
 value: function asec(a) {
 a = Complex.ensureComplex(a);
 var denominator = a[0] ** 2 + a[1] ** 2;

 var realPart = a[0] / denominator;
 var imagPart = -a[1] / denominator;

 return Complex.acos([realPart, imagPart]);
 }
 }, {
 key: "acsc",
 value: function acsc(a) {
 a = Complex.ensureComplex(a);
 var denominator = a[0] ** 2 + a[1] ** 2;

 var realPart = a[0] / denominator;
 var imagPart = -a[1] / denominator;

 return Complex.asin([realPart, imagPart]);
 }
 }, {
 key: "sinh",
 value: function sinh(a) {
 if (Complex.isReal(a)) {
 return Math.sinh(Complex.re(a));
 }
 var realPart = Math.sinh(a[0]) * Math.cos(a[1]);
 var imagPart = Math.cosh(a[0]) * Math.sin(a[1]);
 return [realPart, imagPart];
 }
 }, {
 key: "cosh",
 value: function cosh(a) {
 if (Complex.isReal(a)) {
 return Math.cosh(Complex.re(a));
 }
 var realPart = Math.cosh(a[0]) * Math.cos(a[1]);
 var imagPart = Math.sinh(a[0]) * Math.sin(a[1]);
 return [realPart, imagPart];
 }
 }, {
 key: "tanh",
 value: function tanh(a) {
 if (Complex.isReal(a)) {
 return Math.tanh(Complex.re(a));
 }
 var denominator = Math.cosh(a[0] * 2) + Math.cos(a[1] * 2);
 var realPart = Math.sinh(a[0] * 2) / denominator;
 var imagPart = Math.sin(a[1] * 2) / denominator;
 return [realPart, imagPart];
 }
 }, {
 key: "coth",
 value: function coth(a) {
 a = Complex.ensureComplex(a);
 var denominator = Math.cosh(a[0] * 2) - Math.cos(a[1] * 2);
 var realPart = Math.sinh(a[0] * 2) / denominator;
 var imagPart = -Math.sin(a[1] * 2) / denominator;
 return [realPart, imagPart];
 }
 }, {
 key: "csch",
 value: function csch(a) {
 a = Complex.ensureComplex(a);
 var denominator = Math.cos(a[1] * 2) - Math.cosh(a[0] * 2);
 var realPart = -2 * Math.sinh(a[0]) * Math.cos(a[1]) / denominator;
 var imagPart = 2 * Math.cosh(a[0]) * Math.sin(a[1]) / denominator;
 return [realPart, imagPart];
 }
 }, {
 key: "sech",
 value: function sech(a) {
 a = Complex.ensureComplex(a);
 var denominator = Math.cos(a[1] * 2) + Math.cosh(a[0] * 2);
 var realPart = 2 * Math.cosh(a[0]) * Math.cos(a[1]) / denominator;
 var imagPart = -2 * Math.sinh(a[0]) * Math.sin(a[1]) / denominator;
 return [realPart, imagPart];
 }
 }, {
 key: "asinh",
 value: function asinh(a) {
 if (Complex.isReal(a)) {
 return Math.asinh(Complex.re(a));
 }
 var realPart = void 0,
 imagPart = void 0;

 var _Complex$asin3 = Complex.asin([a[1], -a[0]]);

 var _Complex$asin4 = _slicedToArray(_Complex$asin3, 2);

 realPart = _Complex$asin4[0];
 imagPart = _Complex$asin4[1];

 return [-imagPart, realPart];
 }
 }, {
 key: "acosh",
 value: function acosh(a) {
 if (Complex.isReal(a)) {
 return Math.acosh(Complex.re(a));
 }

 var _Complex$acos = Complex.acos(a),
 _Complex$acos2 = _slicedToArray(_Complex$acos, 2),
 realPart = _Complex$acos2[0],
 imagPart = _Complex$acos2[1];

 if (imagPart <= 0) {
 return [-imagPart, realPart];
 } else {
 return [imagPart, -realPart];
 }
 }
 }, {
 key: "atanh",
 value: function atanh(a) {
 if (Complex.isReal(a)) {
 return Math.atanh(Complex.re(a));
 }
 var denominator = (1 - a[0]) ** 2 + a[1] ** 2;

 var realPart = ((1 + a[0]) * (1 - a[0]) - a[1] ** 2) / denominator;
 var imagPart = (a[1] * (1 - a[0]) + a[1] * (1 + a[0])) / denominator;

 var cache = Math.log(Math.sqrt(realPart ** 2 + imagPart ** 2)) / 2;
 imagPart = Math.atan2(imagPart, realPart) / 2;
 realPart = cache;

 return [realPart, imagPart];
 }
 }, {
 key: "acoth",
 value: function acoth(a) {
 a = Complex.ensureComplex(a);
 var denominator = a[0] ** 2 + a[1] ** 2;

 var realPart = a[0] / denominator;
 var imagPart = -a[1] / denominator;

 return Complex.atanh([realPart, imagPart]);
 }
 }, {
 key: "acsch",
 value: function acsch(a) {
 a = Complex.ensureComplex(a);
 var denominator = a[0] ** 2 + a[1] ** 2;

 var realPart = a[0] / denominator;
 var imagPart = -a[1] / denominator;

 return Complex.asinh([realPart, imagPart]);
 }
 }, {
 key: "asech",
 value: function asech(a) {
 a = Complex.ensureComplex(a);
 var denominator = a[0] ** 2 + a[1] ** 2;

 var realPart = a[0] / denominator;
 var imagPart = -a[1] / denominator;

 return Complex.acosh([realPart, imagPart]);
 }
 }, {
 key: "inverse",
 value: function inverse(a) {
 if (Complex.isReal(a)) {
 return 1 / Complex.re(a);
 }
 a = Complex.ensureComplex(a);
 var denominator = a[0] ** 2 + a[1] ** 2;

 var realPart = a[0] === 0 ? 0 : a[0] / denominator;
 var imagPart = a[1] === 0 ? 0 : -a[1] / denominator;

 return [realPart, imagPart];
 }

 // default imaginary part to zero if non-array is passed

 }, {
 key: "ensureComplex",
 value: function ensureComplex(arg) {
 return arg.constructor === Array ? arg : [arg, 0];
 }

 // Return the real part of complex or just return if float

 }, {
 key: "re",
 value: function re(arg) {
 return arg.constructor === Array ? arg[0] : arg;
 }

 // Return the imaginary part or just 0 if float

 }, {
 key: "im",
 value: function im(arg) {
 return arg.constructor === Array ? arg[1] : 0;
 }
 }, {
 key: "isZero",
 value: function isZero(arg) {
 if (arg.constructor === Array && arg.length > 1) {
 return arg[0] === 0 && arg[1] === 0;
 } else if (arg.constructor === Array && arg.length > 0) {
 return arg[0] === 0;
 } else {
 return arg === 0;
 }
 }
 }, {
 key: "prod",
 value: function prod(numbers) {
 var total = Complex.mul(numbers[0], numbers[1]);
 for (var i = 2; i < numbers.length; i++) {
 total = Complex.mul(total, numbers[i]);
 }
 return total;
 }
 }, {
 key: "sum",
 value: function sum(numbers) {
 var total = [0, 0];
 for (var i = 0; i < numbers.length; i++) {
 if (numbers[i].constructor === Array) {
 total[0] += numbers[i][0];
 total[1] += numbers[i][1];
 } else {
 total[0] += numbers[i];
 }
 }
 return total;
 }

 // Check for size of imaginary component compared to real component.

 }, {
 key: "isReal",
 value: function isReal(x) {
 if (x.constructor === Array) {
 return x[0] + x[1] === x[0];
 } else {
 return true;
 }
 }
 }]);

 return Complex;
}();

exports.default = Complex;
},{}],90:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});

var _createClass = function () { function defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } } return function (Constructor, protoProps, staticProps) { if (protoProps) defineProperties(Constructor.prototype, protoProps); if (staticProps) defineProperties(Constructor, staticProps); return Constructor; }; }();

var _complex = require('../utils/complex.js');

var _complex2 = _interopRequireDefault(_complex);

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

// These can be used for lightweight stuff, but they are much slower (~200 times slower)
// They should compile to c structures well
var ComplexNumber = function () {
 function ComplexNumber() {
 var buffer = arguments.length > 0 && arguments[0] !== undefined ? arguments[0] : new ArrayBuffer(16);

 _classCallCheck(this, ComplexNumber);

 this.buffer = buffer;
 }

 _createClass(ComplexNumber, [{
 key: 'add',


 // example wrapper
 value: function add(complexV) {
 var result = new ComplexNumber();
 _complex2.default.add(this.buffer, 0, complexV.buffer, 0, result.buffer, 0);
 return result;
 }

 // could be created dynamically?

 }, {
 key: 'mul',
 value: function mul(complexV) {
 var result = new ComplexNumber();
 _complex2.default.mul(this.buffer, 0, complexV.buffer, 0, result.buffer, 0);
 return result;
 }

 // about 10 times faster than mul, still 20 x slower than structureless complex numbers

 }, {
 key: 'mulEql',
 value: function mulEql(complexV) {
 _complex2.default.mul(this.buffer, 0, complexV.buffer, 0, this.buffer, 0);
 return this;
 }
 }, {
 key: 'buffer',
 set: function set(buffer) {
 this._buffer = buffer;
 this._r = new Float64Array(this._buffer, 0, 1);
 this._i = new Float64Array(this._buffer, 8, 1);
 },
 get: function get() {
 return this._buffer;
 }
 }, {
 key: 'r',
 set: function set(v) {
 this._r[0] = v;
 },
 get: function get() {
 return this._r[0];
 }
 }, {
 key: 'i',
 set: function set(v) {
 this._i[0] = v;
 },
 get: function get() {
 return this._i[0];
 }
 }]);

 return ComplexNumber;
}();

exports.default = ComplexNumber;
},{"../utils/complex.js":89}],91:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.d1mach = d1mach;
// !***BEGIN PROLOGUE D1MACH
// !***PURPOSE Return floating point machine dependent constants.
// !***LIBRARY SLATEC
// !***CATEGORY R1
// !***TYPE SINGLE PRECISION (D1MACH-S, D1MACH-D)
// !***KEYWORDS MACHINE CONSTANTS
// !***AUTHOR Fox, P. A., (Bell Labs)
// ! Hall, A. D., (Bell Labs)
// ! Schryer, N. L., (Bell Labs)
// !***DESCRIPTION
// !
// ! D1MACH can be used to obtain machine-dependent parameters for the
// ! local machine environment. It is a function subprogram with one
// ! (input) argument, and can be referenced as follows:
// !
// ! A = D1MACH(I)
// !
// ! where I=1,...,5. The (output) value of A above is determined by
// ! the (input) value of I. The results for various values of I are
// ! discussed below.
// !
// ! D1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
// ! D1MACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
// ! D1MACH(3) = B**(-T), the smallest relative spacing.
// ! D1MACH(4) = B**(1-T), the largest relative spacing.
// ! D1MACH(5) = LOG10(B)
// !
// ! Assume single precision numbers are represented in the T-digit,
// ! base-B form
// !
// ! sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
// !
// ! where 0 .LE. X(I) .LT. B for I=1,...,T, 0 .LT. X(1), and
// ! EMIN .LE. E .LE. EMAX.
// !
// ! The values of B, T, EMIN and EMAX are provided in I1MACH as
// ! follows:
// ! I1MACH(10) = B, the base.
// ! I1MACH(11) = T, the number of base-B digits.
// ! I1MACH(12) = EMIN, the smallest exponent E.
// ! I1MACH(13) = EMAX, the largest exponent E.
// !
// !
// !***REFERENCES P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
// ! a portable library, ACM Transactions on Mathematical
// ! Software 4, 2 (June 1978), pp. 177-188.
// !***ROUTINES CALLED XERMSG
// !***REVISION HISTORY (YYMMDD)
// ! 790101 DATE WRITTEN
// ! 960329 Modified for Fortran 90 (BE after suggestions by EHG)
// !***END PROLOGUE D1MACH
// !
function d1mach(input) {
 switch (input) {
 case 1:
 return Number.MIN_VALUE;
 case 2:
 return Number.MAX_VALUE;
 case 3:
 return Number.EPSILON / 2; // the smallest relative spacing.
 case 4:
 return Number.EPSILON; // the largest relative spacing.
 case 5:
 return Math.log10(2);
 default:
 throw new Error('d1mach expects an integer 1-5.');
 }
}
},{}],92:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
exports.i1mach = i1mach;
// import * as py from '../../utils/pythonHelpers.js';
// !DECK I1MACH
// INTEGER FUNCTION I1MACH (I)
// IMPLICIT NONE
// INTEGER :: I
// REAL :: X
// DOUBLE PRECISION :: XX
// !***BEGIN PROLOGUE I1MACH
// !***PURPOSE Return integer machine dependent constants.
// !***LIBRARY SLATEC
// !***CATEGORY R1
// !***TYPE INTEGER (I1MACH-I)
// !***KEYWORDS MACHINE CONSTANTS
// !***AUTHOR Fox, P. A., (Bell Labs)
// ! Hall, A. D., (Bell Labs)
// ! Schryer, N. L., (Bell Labs)
// !***DESCRIPTION
// !
// ! I1MACH can be used to obtain machine-dependent parameters for the
// ! local machine environment. It is a function subprogram with one
// ! (input) argument and can be referenced as follows:
// !
// ! K = I1MACH(I)
// !
// ! where I=1,...,16. The (output) value of K above is determined by
// ! the (input) value of I. The results for various values of I are
// ! discussed below.
// !
// ! I/O unit numbers:
// ! I1MACH( 1) = the standard input unit.
// ! I1MACH( 2) = the standard output unit.
// ! I1MACH( 3) = the standard punch unit.
// ! I1MACH( 4) = the standard error message unit.
// !
// ! Words:
// ! I1MACH( 5) = the number of bits per integer storage unit.
// ! I1MACH( 6) = the number of characters per integer storage unit.
// !
// ! Integers:
// ! assume integers are represented in the S-digit, base-A form
// !
// ! sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
// !
// ! where 0 .LE. X(I) .LT. A for I=0,...,S-1.
// ! I1MACH( 7) = A, the base.
// ! I1MACH( 8) = S, the number of base-A digits.
// ! I1MACH( 9) = A**S - 1, the largest magnitude.
// !
// ! Floating-Point Numbers:
// ! Assume floating-point numbers are represented in the T-digit,
// ! base-B form
// ! sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
// !
// ! where 0 .LE. X(I) .LT. B for I=1,...,T,
// ! 0 .LT. X(1), and EMIN .LE. E .LE. EMAX.
// ! I1MACH(10) = B, the base.
// !
// ! Single-Precision:
// ! I1MACH(11) = T, the number of base-B digits.
// ! I1MACH(12) = EMIN, the smallest exponent E.
// ! I1MACH(13) = EMAX, the largest exponent E.
// !
// ! Double-Precision:
// ! I1MACH(14) = T, the number of base-B digits.
// ! I1MACH(15) = EMIN, the smallest exponent E.
// ! I1MACH(16) = EMAX, the largest exponent E.
// !
// ! To alter this function for a particular environment, the desired
// ! set of DATA statements should be activated by removing the C from
// ! column 1. Also, the values of I1MACH(1) - I1MACH(4) should be
// ! checked for consistency with the local operating system.
// !
// !***REFERENCES P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
// ! a portable library, ACM Transactions on Mathematical
// ! Software 4, 2 (June 1978), pp. 177-188.
// !***ROUTINES CALLED (NONE)
// !***REVISION HISTORY (YYMMDD)
// ! 750101 DATE WRITTEN
// ! 960411 Modified for Fortran 90 (BE after suggestions by EHG).
// ! 980727 Modified value of I1MACH(6) (BE after suggestion by EHG).
// !***END PROLOGUE I1MACH
// !
function i1mach(input) {
 var bitSize = 64; // JS is always 64bit
 // The following can be computed with the help of the frexp function in
 // utils/pythonHelpers:
 var digits = 53; // py.frexp(Number.MAX_SAFE_INTEGER)[1];
 var maxExponent = 1024; // py.frexp(Number.MAX_VALUE)[1];
 var minExponent = -1073; // py.frexp(Number.MIN_VALUE)[1];
 switch (input) {
 case 1:
 return 5;
 case 2:
 return 6;
 case 3:
 return 0;
 case 4:
 return 0;
 case 5:
 return bitSize;
 case 6:
 return 4;
 case 7:
 return 2;
 case 8:
 return bitSize - 1;
 case 9:
 return Number.MAX_VALUE;
 case 10:
 return 2;
 case 11:
 return digits;
 case 12:
 return minExponent;
 case 13:
 return maxExponent;
 case 14:
 return digits;
 case 15:
 return minExponent;
 case 16:
 return maxExponent;
 default:
 throw new Error('i1mach expects an integer 1-16.');
 }
}
},{}],93:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
// https://gcc.gnu.org/onlinedocs/gcc-3.4.4/g77/Sign-Intrinsic.html#Sign-Intrinsic
function sign(a, b) {
 var s = b >= 0 ? 1 : -1;
 return Math.abs(a) * s;
}

exports.sign = sign;
},{}],94:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});
// Python has some built-in stuff that javascript doesn't.

// frexp and ldexp from:
// http://croquetweak.blogspot.com/2014/08/deconstructing-floats-frexp-and-ldexp.html
function frexp(value) {
 if (value === 0) return [value, 0];
 var data = new DataView(new ArrayBuffer(8));
 data.setFloat64(0, value);
 var bits = data.getUint32(0) >>> 20 & 0x7FF;
 if (bits === 0) {
 // denormal
 data.setFloat64(0, value * Math.pow(2, 64)); // exp + 64
 bits = (data.getUint32(0) >>> 20 & 0x7FF) - 64;
 }
 var exponent = bits - 1022;
 var mantissa = ldexp(value, -exponent);
 return [mantissa, exponent];
}

function ldexp(mantissa, exponent) {
 var steps = Math.min(3, Math.ceil(Math.abs(exponent) / 1023));
 var result = mantissa;
 for (var i = 0; i < steps; i++) {
 result *= Math.pow(2, Math.floor((exponent + i) / steps));
 }
 return result;
}

function range(start, end, step) {
 var _end = end || start;
 var _start = end ? start : 0;
 var _step = step || 1;
 return Array((_end - _start) / _step).fill(0).map(function (v, i) {
 return _start + i * _step;
 });
}

// * It doesn't handle variadic arguments, that could be supported by js has some perf
// issues with ... (so I've heard). That could be handled like
// zip([[arr1], [arr2], [arr3]]) and then:
// (...rows) => [...rows[0]].map((_,c) => rows.map(row => row[c]))
// * It isn't it's own inverse like the python version. For that, use unzip
function zip(arr1, arr2) {
 return arr1.map(function (_, index) {
 return [arr1[index], arr2[index]];
 });
}

// Does the reverse of zip: [letters, numbers] = unzip( [['a',1], ['b', 2]] )
// letters = ['a', 'b']
// numbers = [1, 2]
function unzip(arr) {
 var res0 = [];
 var res1 = [];
 arr.map(function (pair) {
 res0.push(pair[0]);
 res1.push(pair[1]);
 });
 return [res0, res1];
}

// It is common to want to combine zip with the python pattern that iterates and creates arrays ([func(a) for a in alpha_s] similar to `map` in js). So this function combines the two a little more smoothly.
function zipmap(arr1, arr2, func) {
 return arr1.map(function (_, index) {
 return func(arr1[index], arr2[index]);
 });
}

// sum all elements of an array. This is much faster than reduce: https://jsperf.com/js-sum-3367890876
function sum(numbers) {
 var total = 0;
 for (var i = 0; i < numbers.length; i++) {
 total += numbers[i];
 }
 return total;
}

// Mimick python's get method for easier translation
function get(obj, key) {
 var defaultValue = arguments.length > 2 && arguments[2] !== undefined ? arguments[2] : null;

 if (obj.hasOwnProperty(key)) {
 return obj[key];
 } else {
 return defaultValue;
 }
}

function mod(x, y) {
 return (x % y + y) % y;
}

// divmod
function divmod(x, y) {
 return [Math.floor(x / y), mod(x, y)];
}

exports.frexp = frexp;
exports.ldexp = ldexp;
exports.range = range;
exports.zip = zip;
exports.zipmap = zipmap;
exports.unzip = unzip;
exports.sum = sum;
exports.get = get;
exports.mod = mod;
exports.divmod = divmod;
},{}],95:[function(require,module,exports){
Object.defineProperty(exports, "__esModule", {
 value: true
});

var _chebyshev = require('./lib/special/chebyshev.js');

var _chebyshev2 = _interopRequireDefault(_chebyshev);

var _complex = require('./lib/utils/complex.js');

var _complex2 = _interopRequireDefault(_complex);

var _complexNumber = require('./lib/utils/complexNumber.js');

var _complexNumber2 = _interopRequireDefault(_complexNumber);

var _cephes = require('./lib/special/cephes.js');

var _cephes2 = _interopRequireDefault(_cephes);

var _bessel = require('./lib/special/mpmath/functions/bessel.js');

var _quadpack = require('./lib/integration/quadpack.js');

var _amos = require('./lib/special/amos.js');

var amos = _interopRequireWildcard(_amos);

var _struve = require('./lib/special/c_misc/struve.js');

var struve = _interopRequireWildcard(_struve);

function _interopRequireWildcard(obj) { if (obj && obj.__esModule) { return obj; } else { var newObj = {}; if (obj != null) { for (var key in obj) { if (Object.prototype.hasOwnProperty.call(obj, key)) newObj[key] = obj[key]; } } newObj.default = obj; return newObj; } }

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

exports.default = {
 Chebyshev: _chebyshev2.default,
 Complex: _complex2.default,
 ComplexNumber: _complexNumber2.default,
 Cephes: _cephes2.default,
 besseli: _bessel.besseli,
 besselk: _bessel.besselk,
 struvel: _bessel.struvel,
 quad: _quadpack.quad,
 amos: amos,
 struve: struve
};
},{"./lib/integration/quadpack.js":2,"./lib/special/amos.js":9,"./lib/special/c_misc/struve.js":47,"./lib/special/cephes.js":48,"./lib/special/chebyshev.js":82,"./lib/special/mpmath/functions/bessel.js":85,"./lib/utils/complex.js":89,"./lib/utils/complexNumber.js":90}]},{},[1]);
;
 let element = args[args.length - 3]
 let step = args[args.length - 2]
 let n = args[args.length - 1]
 let results = []
 // console.log('emcoil:', args[0], `element: ${element}, step: ${step}, n: ${n}`);
 progress(0);
 for (let i = 0; i < n; i++) {
   progress(i/n)
   let fnArgs = [...args[0]]
   fnArgs[element] += step * i
   // console.log(i, '_fnArgs:', fnArgs)
   let result = [ fnArgs[0], fnArgs[1], (axialForce(...fnArgs) || 0), (radialForce(...fnArgs) || 0) ]
   results.push(result)
 }
 progress(1)
 return results
 }
