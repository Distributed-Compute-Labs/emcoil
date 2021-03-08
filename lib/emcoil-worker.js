function workFn(...args) {

  // Import from SparcMath
  let SparcMath = require('dcp-math');

  // console.log(SparcMath);
  let besselI = SparcMath.default.besseli;
  let struveL = SparcMath.default.struvel;
  let Complex = SparcMath.default.Complex;
  let quad = SparcMath.default.quad;
  let qag = SparcMath.default.qag;
  let amos = SparcMath.default.amos;
  let besselK = amos.besselK;

  // Note to dan you can ignore these constants
  // just added them so I could search and replace them
  // to get that table in the current demo working
  // just make whatever you want

  // const _$_a = 0;
  // const _$_b = 0;
  // const _$_l = 0;
  // const _$_n = 0;
  // const _$_i = 0;
  // const _$_v = 0;
  // const _$_R_i = 0;
  // const _$_R_o = 0;
  // const _$_sigma = 0;

  // Simulation parameters
  // const a = _$_a || 0.01; // m - coil inner radius
  // const b = _$_b || 0.02; // m - coil outer radius
  // const length = _$_l || 0.01; // m - coil length
  // const n = _$_n || 1e5; // m^-2 - turn density
  // const i = _$_i || 1; // A - coil current
  // const nu = _$_v || 10; // m * s^-1 - coil velocity
  // const rInner = _$_R_i || 0.03; // m - tube inner radius
  // const rOuter = _$_R_o || 0.04; // m - tube outer radius
  // const sigma = _$_sigma || 5.96e7; // S * m^-1 - tube conductivity

  const mu0 = 1.25664e-6; // permeability


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


  {
    // Note to Dan, this is another variable that I search and replace
    var $$body;

    let element = args[args.length - 3]
    let step = args[args.length - 2]
    let num = args[args.length - 1]

    let results = []

    console.log('emcoil:', args[0], `element: ${element}, step: ${step}, num: ${num}`);
    progress(0);
    for (let i = 0; i < num; i++) {
      progress(i/num)
      let fnArgs = [...args[0]]
      fnArgs[element] += step * i
      // console.log(i, '_fnArgs:', fnArgs)
      let result = [ fnArgs[0], fnArgs[1], (axialForce(...fnArgs) || 0), (radialForce(...fnArgs) || 0) ]
      results.push(result)
    }
    progress(1)
    console.log('emcoil: returning results', results.length);
    return results
  }
}
