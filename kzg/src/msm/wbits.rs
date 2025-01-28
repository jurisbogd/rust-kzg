use crate::{Fr, G1Affine, G1Fp, G1GetFp, G1Mul, Scalar256, G1};

use blst::{
    blst_p1, blst_p1_affine, blst_p1s_mult_wbits_precompute_sizeof, blst_p1s_to_affine, byte,
};
use core::{marker::PhantomData, ptr::null, slice};

use super::pippenger_utils;

const WBITS: usize = 2;

#[derive(Debug, Clone)]
pub struct WbitsTable<
    TFr: Fr,
    TG1: G1 + G1Mul<TFr> + G1GetFp<TG1Fp>,
    TG1Fp: G1Fp,
    TG1Affine: G1Affine<TG1, TG1Fp>,
> {
    precomputation_table: Vec<TG1Affine>,

    g1_marker: PhantomData<TG1>,
    g1_fp_marker: PhantomData<TG1Fp>,
    fr_marker: PhantomData<TFr>,
    g1_affine_marker: PhantomData<TG1Affine>,
}

impl<
        TFr: Fr,
        TG1Fp: G1Fp,
        TG1: G1 + G1Mul<TFr> + G1GetFp<TG1Fp>,
        TG1Affine: G1Affine<TG1, TG1Fp>,
    > WbitsTable<TFr, TG1, TG1Fp, TG1Affine>
{
    pub fn new(points: &[TG1]) -> Result<Option<Self>, String> {
        //Allocate memory for precomputation table
        let precomputation_table_size = Self::get_precomputation_table_size(points.len());
        let mut precomputation_table = vec![TG1Affine::default(); precomputation_table_size];

        unsafe {
            //Convert G1 points to blst_p1_affine points
            let points_arr = [points.as_ptr() as *const blst_p1, null()];
            let mut blst_p1_affines = vec![blst_p1_affine::default(); points.len()];
            blst_p1s_to_affine(
                blst_p1_affines.as_mut_ptr(),
                points_arr.as_ptr(),
                points.len(),
            );

            g1s_precompute_wbits(
                precomputation_table.as_mut_ptr() as *mut TG1Affine,
                WBITS,
                blst_p1_affines.as_ptr() as *const TG1Affine,
                points.len(),
            );
        }

        Ok(Some(Self {
            precomputation_table,

            fr_marker: PhantomData,
            g1_marker: PhantomData,
            g1_fp_marker: PhantomData,
            g1_affine_marker: PhantomData,
        }))
    }

    pub fn multiply_sequential(&self, tfrs: &[TFr]) -> TG1 {
        let mut ret = TG1::default();
        let scalars = tfrs
            .into_iter()
            .map(|tfr| tfr.to_scalar())
            .collect::<Vec<_>>();
        let mut scratch = vec![TG1::default(); 8192];

        unsafe {
            mult_wbits_simple_scalars(
                &mut ret as *mut TG1,
                self.precomputation_table.as_ptr() as *const TG1Affine,
                WBITS,
                scalars.len(),
                scalars.as_ptr() as *const byte,
                size_of::<TFr>() * 8,
                scratch.as_mut_ptr() as *mut TG1,
            );
        }

        ret
    }

    #[cfg(feature = "parallel")]
    pub fn multiply_parallel(&self, scalars: &[TFr]) -> TG1 {
        self.multiply_sequential(scalars)
    }

    fn get_precomputation_table_size(npoints: usize) -> usize {
        npoints << (WBITS - 1)
    }
}

unsafe fn g1s_precompute_wbits<
    TG1: G1 + G1GetFp<TG1Fp>,
    TG1Fp: G1Fp,
    TG1Affine: G1Affine<TG1, TG1Fp>,
>(
    table: *mut TG1Affine,
    wbits: usize,
    mut points: *const TG1Affine,
    mut npoints: usize,
) {
    //size_t total = npoints << (wbits-1);
    let total = npoints << (wbits - 1);
    //size_t nwin = (size_t)1 << (wbits-1);
    let nwin = (1 as usize) << (wbits - 1);
    //size_t nmin = wbits>9 ? (size_t)1: (size_t)1 << (9-wbits);
    let nmin = if wbits > 9 {
        1 as usize
    } else {
        (1 as usize) << (9 - wbits)
    };
    //size_t i, top = 0;
    let mut i;
    let mut top = 0;
    //ptype *rows, *row;
    let mut rows;
    let mut row;
    //const ptype##_affine *point = NULL;
    let mut point;
    //size_t stride = ((512*1024)/sizeof(ptype##_affine)) >> wbits;
    let mut stride = (((512 as usize) * (1024 as usize)) / size_of::<blst_p1_affine>()) >> wbits;
    //if (stride == 0) stride = 1;
    if stride == 0 {
        stride = 1;
    }

    //while (npoints >= nmin) {
    'outer: while npoints >= nmin {
        //size_t limit = total - npoints;
        let limit = total - npoints;

        //if (top + (stride << wbits) > limit) {
        if (top + (stride << wbits)) > limit {
            //stride = (limit - top) >> wbits;
            stride = (limit - top) >> wbits;
            //if (stride == 0) break;
            if stride == 0 {
                break 'outer;
            }
        }
        //rows = row = (ptype *)(&table[top]);
        row = (table.add(top)) as *mut TG1;
        rows = row;
        //for (i = 0; i < stride; i++, row += nwin) \
        for _ in 0..stride {
            //point = *points ? *points++ : point+1,
            point = points;
            points = points.add(1);
            let row_slice = slice::from_raw_parts_mut(row, nwin);
            //ptype##_precompute_row_wbits(row, wbits, point);
            g1_precompute_row_wbits(row_slice, &*point);

            row = row.add(nwin);
        }
        // g1s_to_affine_row_wbits(&table[top], rows, wbits, npoints);
        g1s_to_affine_row_wbits(table.add(top), rows, wbits, stride);
        //top += stride << (wbits-1);
        top += stride << (wbits - 1);
        // npoints -= stride;
        npoints -= stride;
    }
    //rows = row = alloca(2*sizeof(ptype##_affine) * npoints * nwin);
    let mut alloca = vec![0 as byte; 2 * size_of::<blst_p1_affine>() * npoints * nwin];
    //let mut alloca = vec![TG1Affine::default(); 2 * npoints * nwin];
    row = alloca.as_mut_ptr() as *mut TG1;
    rows = row;
    //for (i = 0; i < npoints; i++, row += nwin)
    i = 0;
    while i < npoints {
        //point = *points ? *points++ : point+1,
        point = points;
        points = points.add(1);
        let row_slice = slice::from_raw_parts_mut(row, nwin);
        //ptype##_precompute_row_wbits(row, wbits, point);
        g1_precompute_row_wbits(row_slice, &*point);

        i += 1;
        row = row.add(nwin);
    }
    //ptype##s_to_affine_row_wbits(&table[top], rows, wbits, npoints);
    g1s_to_affine_row_wbits(table.add(top), rows, wbits, npoints);
}

fn g1_precompute_row_wbits<
    TG1: G1 + G1GetFp<TG1Fp>,
    TG1Fp: G1Fp,
    TG1Affine: G1Affine<TG1, TG1Fp>,
>(
    row: &mut [TG1],
    point: &TG1Affine,
) {
    *row[0].x_mut() = *point.x();
    *row[0].y_mut() = *point.y();
    *row[0].z_mut() = G1Fp::one();

    // row[1]=p*(1+1)
    row[1] = row[0].dbl();

    for i in (2..row.len()).step_by(2) {
        // row[2]=p*(2+1)
        row[i] = add_affine(&row[i - 1], point);
        // row[3]=p*(2+2)
        row[i + 1] = row[i / 2].dbl();
    }
}

// size_t prefix##s_mult_wbits_precompute_sizeof(size_t wbits, size_t npoints) \
// { return (sizeof(ptype##_affine)*npoints) << (wbits-1); } \
// fn get_precomputation_table_size<TG1: G1, TG1Fp: G1Fp, TG1Affine: G1Affine<TG1, TG1Fp>>(
//     npoints: usize,
// ) -> usize {
//     size_of::<TG1Affine>() * npoints << (WBITS - 1)
// }

fn add_affine<TG1: G1 + G1GetFp<TG1Fp>, TG1Fp: G1Fp, TG1Affine: G1Affine<TG1, TG1Fp>>(
    point: &TG1,
    affine: &TG1Affine,
) -> TG1 {
    //ptype p3;
    let mut p3 = TG1::default();
    //vec##bits Z1Z1, H, HH, I, J;
    let z1z1;
    let mut h;
    let hh;
    let mut i;
    let j;

    if point.z().is_zero() {
        let mut out = TG1::default();
        *out.x_mut() = *affine.x();
        *out.y_mut() = *affine.y();
        *out.z_mut() = TG1Fp::one();
        return out;
    }
    if affine.x().is_zero() && affine.y().is_zero() {
        return point.clone();
    }

    // sqr_##field(Z1Z1, p1->Z);           /* Z1Z1 = Z1^2 */
    z1z1 = point.z().square();

    // mul_##field(p3.Z, Z1Z1, p1->Z);     /* Z1*Z1Z1 */
    *p3.z_mut() = z1z1.mul_fp(point.z());
    // mul_##field(p3.Z, p3.Z, affine->Y);     /* S2 = Y2*Z1*Z1Z1 */
    *p3.z_mut() = p3.z().mul_fp(affine.y());

    // mul_##field(H, affine->X, Z1Z1);        /* U2 = X2*Z1Z1 */
    h = affine.x().mul_fp(&z1z1);
    // sub_##field(H, H, p1->X);           /* H = U2-X1 */
    h = h.sub_fp(point.x());

    // sqr_##field(HH, H);                 /* HH = H^2 */
    hh = h.square();
    // add_##field(I, HH, HH);
    i = hh.add_fp(&hh);
    // add_##field(I, I, I);               /* I = 4*HH */
    i = i.add_fp(&i);

    // mul_##field(p3.Y, p1->X, I);        /* V = X1*I */
    *p3.y_mut() = point.x().mul_fp(&i);
    // mul_##field(J, H, I);               /* J = H*I */
    j = h.mul_fp(&i);
    // mul_##field(I, J, p1->Y);           /* Y1*J */
    i = j.mul_fp(point.y());

    // sub_##field(p3.Z, p3.Z, p1->Y);     /* S2-Y1 */
    *p3.z_mut() = p3.z().sub_fp(point.y());
    // add_##field(p3.Z, p3.Z, p3.Z);      /* r = 2*(S2-Y1) */
    *p3.z_mut() = p3.z().add_fp(p3.z());

    // sqr_##field(p3.X, p3.Z);            /* r^2 */
    *p3.x_mut() = p3.z().square();
    // sub_##field(p3.X, p3.X, J);         /* r^2-J */
    *p3.x_mut() = p3.x().sub_fp(&j);
    // sub_##field(p3.X, p3.X, p3.Y);
    *p3.x_mut() = p3.x().sub_fp(p3.y());
    // sub_##field(p3.X, p3.X, p3.Y);      /* X3 = r^2-J-2*V */
    *p3.x_mut() = p3.x().sub_fp(p3.y());

    // sub_##field(p3.Y, p3.Y, p3.X);      /* V-X3 */
    *p3.y_mut() = p3.y().sub_fp(p3.x());
    // mul_##field(p3.Y, p3.Y, p3.Z);      /* r*(V-X3) */
    *p3.y_mut() = p3.y().mul_fp(p3.z());
    // sub_##field(p3.Y, p3.Y, I);
    *p3.y_mut() = p3.y().sub_fp(&i);
    // sub_##field(p3.Y, p3.Y, I);         /* Y3 = r*(V-X3)-2*Y1*J */
    *p3.y_mut() = p3.y().sub_fp(&i);

    // add_##field(p3.Z, p1->Z, H);        /* Z1+H */
    *p3.z_mut() = point.z().add_fp(&h);
    // sqr_##field(p3.Z, p3.Z);            /* (Z1+H)^2 */
    *p3.z_mut() = p3.z().square();
    // sub_##field(p3.Z, p3.Z, Z1Z1);      /* (Z1+H)^2-Z1Z1 */
    *p3.z_mut() = p3.z().sub_fp(&z1z1);
    // sub_##field(p3.Z, p3.Z, HH);        /* Z3 = (Z1+H)^2-Z1Z1-HH */
    *p3.z_mut() = p3.z().sub_fp(&hh);

    return p3;
}

unsafe fn g1s_to_affine_row_wbits<
    TG1: G1 + G1GetFp<TG1Fp>,
    TG1Fp: G1Fp,
    TG1Affine: G1Affine<TG1, TG1Fp>,
>(
    dst: *mut TG1Affine,
    src: *mut TG1,
    wbits: usize,
    npoints: usize,
) {
    // POINTonE1s_to_affine_row_wbits(
    //     dst as *mut blst_p1_affine,
    //     src as *mut blst_p1,
    //     wbits,
    //     npoints,
    // );
    //to_affine_row_wbits_no_ptrs(dst, src, wbits, npoints);
    to_affine_row_wbits(dst, src, wbits, npoints);
}

unsafe fn to_affine_row_wbits<
    TG1: G1 + G1GetFp<TG1Fp>,
    TG1Fp: G1Fp,
    TG1Affine: G1Affine<TG1, TG1Fp>,
>(
    mut dst: *mut TG1Affine,
    mut src: *mut TG1,
    wbits: usize,
    npoints: usize,
) {
    let total = npoints << (wbits - 1);
    let nwin = (1 as usize) << (wbits - 1);
    let mut acc;
    let mut zz;
    let mut zzz;

    //src += total;
    src = src.add(total);
    //acc = (vec##bits *)src;
    acc = src as *mut TG1Fp;
    //vec_copy(acc++, one, sizeof(vec##bits));
    *acc = TG1Fp::one();
    acc = acc.add(1);

    {
        // //for (i = 0; i < npoints; i++)
        for _ in 0..npoints {
            //     //for (j = nwin; --src, --j; acc++)
            for _ in 1..nwin {
                src = src.sub(1);

                *acc = (*acc.sub(1)).mul_fp((*src).z());

                acc = acc.add(1);
            }
            src = src.sub(1);
        }
    }

    //--acc; reciprocal_##field(acc[0], acc[0]);
    acc = acc.sub(1);
    *acc = (*acc).inverse().unwrap();
    //blst_fp_inverse(acc as *mut blst_fp, acc as *const blst_fp);

    {
        let mut i = 0;
        while i < npoints {
            //vec_copy(dst++, src++, sizeof(ptype##_affine));
            *(*dst).x_mut() = *(*src).x();
            *(*dst).y_mut() = *(*src).y();
            dst = dst.add(1);
            src = src.add(1);
            //for (j = 1; j < nwin; j++, acc--, src++, dst++) {
            let mut j = 1;
            while j < nwin {
                let acc_s1 = acc.sub(1);

                //mul_##field(acc[-1], acc[-1], acc[0]);  /* 1/Z        */
                //let mut bzz = (*acc.sub(1)).sub_fp(&*acc);
                *acc_s1 = (*acc_s1).mul_fp(&*acc);
                //sqr_##field(ZZ, acc[-1]);               /* 1/Z^2      */
                zz = (*acc_s1).square();
                //mul_##field(ZZZ, ZZ, acc[-1]);          /* 1/Z^3      */
                zzz = zz.mul_fp(&*acc_s1);
                //mul_##field(acc[-1], src->Z, acc[0]);
                *acc_s1 = (*src).z().mul_fp(&*acc);
                // mul_##field(dst->X, src->X, ZZ);        /* X = X'/Z^2 */
                *(*dst).x_mut() = (*src).x().mul_fp(&zz);
                // mul_##field(dst->Y, src->Y, ZZZ);       /* Y = Y'/Z^3 */
                *(*dst).y_mut() = (*src).y().mul_fp(&zzz);

                j += 1;
                acc = acc.sub(1);
                src = src.add(1);
                dst = dst.add(1);
            }

            i += 1;
        }
    }
}

type LimbT = u64;

const SCRATCH_SZ: usize = 8192;

// void ptype##s_mult_wbits_simple_scalars(ptype *ret, const ptype##_affine table[], \
//                                         size_t wbits, size_t npoints, \
// 										const byte *scalars, size_t nbits, \
// 										ptype scratch[]) \
unsafe fn mult_wbits_simple_scalars<
    TG1: G1 + G1GetFp<TG1Fp>,
    TG1Fp: G1Fp,
    TG1Affine: G1Affine<TG1, TG1Fp>,
>(
    ret: *mut TG1,
    table: *const TG1Affine,
    wbits: usize,
    npoints: usize,
    scalars: *const byte,
    mut nbits: usize,
    scratch: *mut TG1,
) {
    // limb_t wmask, wval; \
    let mut wmask;
    // size_t i, j, z, nbytes, window, nwin = (size_t)1 << (wbits-1); \
    let mut i;
    let mut j;
    let nwin = (1 as usize) << (wbits - 1);
    // const byte *scalar = scalars; \
    let mut scalar = scalars;
    // const ptype##_affine *row = table; \
    let mut row = table;

    // nbytes = (nbits + 7)/8; /* convert |nbits| to bytes */ \
    let nbytes = (nbits + 7) / 8;

    // /* top excess bits modulo target window size */ \
    // window = nbits % wbits; /* yes, it may be zero */ \
    let mut window = nbits % wbits;
    // wmask = ((limb_t)1 << (window + 1)) - 1; \
    wmask = ((1 as LimbT) << (window + 1)) - 1;

    // nbits -= window; \
    nbits -= window;
    // z = is_zero(nbits); \
    let z = if nbits == 0 { 1 } else { 0 };

    {
        // wval = (get_wval_limb(scalar, nbits - (z^1), wbits + (z^1)) << z) & wmask; \
        let mut wval =
            (get_wval_limb_wrapper(scalar, nbits - (z ^ 1), wbits + (z ^ 1)) << z) & wmask;
        // wval = booth_encode(wval, wbits); \
        wval = booth_encode(wval, wbits);
        // ptype##_gather_booth_wbits(&scratch[0], row, wbits, wval); \
        gather_booth_wbits_wrapper(scratch, row, wbits, wval);
        // row += nwin; \
        row = row.add(nwin);
        // scalar += nbytes; \
        scalar = scalar.add(nbytes);
    }

    // i = 1; vec_zero(ret, sizeof(*ret)); \
    i = 1;
    *ret = TG1::zero();
    // while (nbits > 0) { \
    while nbits > 0 {
        // for (j = i; i < npoints; i++, j++, row += nwin, scalar += nbytes) { \
        j = i;
        while i < npoints {
            // if (j == SCRATCH_SZ(ptype)) { \
            if j == SCRATCH_SZ {
                // ptype##s_accumulate(ret, scratch, j); \
                accumulate_wrapper::<TG1, TG1Fp, TG1Affine>(ret, scratch, j);
                // j = 0; \
                j = 0;
            }

            // wval = get_wval_limb(scalar, nbits - 1, window + 1) & wmask; \
            let mut wval = get_wval_limb_wrapper(scalar, nbits - 1, window + 1) & wmask;
            // wval = booth_encode(wval, wbits); \
            wval = booth_encode(wval, wbits);
            // ptype##_gather_booth_wbits(&scratch[j], row, wbits, wval); \
            gather_booth_wbits_wrapper(scratch.add(j), row, wbits, wval);

            i += 1;
            j += 1;
            row = row.add(nwin);
            scalar = scalar.add(nbytes);
        }

        // ptype##s_accumulate(ret, scratch, j); \
        accumulate_wrapper::<TG1, TG1Fp, TG1Affine>(ret, scratch, j);

        // for (j = 0; j < wbits; j++) \
        j = 0;
        while j < wbits {
            // ptype##_double(ret, ret); \
            *ret = (*ret).dbl();

            j += 1;
        }

        // window = wbits; \
        window = wbits;
        // wmask = ((limb_t)1 << (window + 1)) - 1; \
        wmask = ((1 as LimbT) << (window + 1)) - 1;
        // nbits -= window; \
        nbits -= window;
        // i = 0; row = table; \
        i = 0;
        row = table;
        // scalar = scalars; \
        scalar = scalars;
    }

    // for (j = i; i < npoints; i++, j++, row += nwin, scalar += nbytes) { \
    j = i;
    while i < npoints {
        // if (j == SCRATCH_SZ(ptype)) { \
        if j == SCRATCH_SZ {
            // ptype##s_accumulate(ret, scratch, j); \
            accumulate_wrapper::<TG1, TG1Fp, TG1Affine>(ret, scratch, j);
            // j = 0; \
            j = 0;
        }

        // wval = (get_wval_limb(scalar, 0, wbits) << 1) & wmask; \
        let mut wval = (get_wval_limb_wrapper(scalar, 0, wbits) << 1) & wmask;
        // wval = booth_encode(wval, wbits); \
        wval = booth_encode(wval, wbits);
        // ptype##_gather_booth_wbits(&scratch[j], row, wbits, wval); \
        gather_booth_wbits_wrapper(scratch.add(j), row, wbits, wval);

        i += 1;
        j += 1;
        row = row.add(nwin);
        scalar = scalar.add(nbytes);
    }
    // ptype##s_accumulate(ret, scratch, j); \
    accumulate_wrapper::<TG1, TG1Fp, TG1Affine>(ret, scratch, j);
}

fn gather_booth_wbits_wrapper<
    TG1: G1 + G1GetFp<TG1Fp>,
    TG1Fp: G1Fp,
    TG1Affine: G1Affine<TG1, TG1Fp>,
>(
    p: *mut TG1,
    row: *const TG1Affine,
    wbits: usize,
    booth_idx: LimbT,
) {
    unsafe {
        let point = &mut *p;
        let nwin = (1 as usize) << (wbits - 1);
        let row_slice = slice::from_raw_parts(row, nwin);
        gather_booth_wbits(point, row_slice, wbits, booth_idx);
    }
}

unsafe fn gather_booth_wbits<
    TG1: G1 + G1GetFp<TG1Fp>,
    TG1Fp: G1Fp,
    TG1Affine: G1Affine<TG1, TG1Fp>,
>(
    point: &mut TG1,
    row: &[TG1Affine],
    wbits: usize,
    mut booth_idx: LimbT,
) {
    // bool_t booth_sign = (booth_idx >> wbits) & 1; \
    let booth_sign = (booth_idx >> wbits) & 1;
    // static const ptype##_affine infinity = { 0 }; \
    let infinity = TG1Affine::zero();

    // booth_idx &= ((limb_t)1 << wbits) - 1; \
    booth_idx &= ((1 as LimbT) << wbits) - 1;
    // idx_is_zero = is_zero(booth_idx); \
    let idx_is_zero = if booth_idx == 0 { 1 } else { 0 };
    // booth_idx -= 1 ^ idx_is_zero; \
    booth_idx -= 1 ^ idx_is_zero;
    // vec_select(p, &infinity, &row[booth_idx], sizeof(row[0]), idx_is_zero); \
    if idx_is_zero == 1 {
        *point.x_mut() = *infinity.x();
        *point.y_mut() = *infinity.y();
    } else {
        let booth_point = row[booth_idx as usize];
        *point.x_mut() = *booth_point.x();
        *point.y_mut() = *booth_point.y();
    }
    // ptype##_cneg(p, booth_sign); \
    if booth_sign == 1 {
        //cneg_fp(p->Y, p->Y, cbit);
        point.y_mut().neg_assign();
    }
}

fn get_wval_limb_wrapper(d: *const u8, off: usize, bits: usize) -> LimbT {
    if true {
        unsafe { get_wval_limb(d, off, bits) }
    } else {
        let scalar = unsafe { &*(d as *const Scalar256) };
        //get_wval_limb_mine(scalar, off, bits)
        pippenger_utils::get_wval_limb(scalar, off, bits)
    }
}

unsafe fn get_wval_limb(d: *const u8, offset: usize, bits: usize) -> LimbT {
    let mut bytes = d.add(offset / 8);

    let top = (offset + bits - 1) / 8 - offset / 8 + 1;

    let mut out = *bytes as LimbT;
    for i in 1..4 {
        if i < top {
            bytes = bytes.add(1);
            out |= (*bytes as LimbT) << (i * 8);
        }
    }

    out >> (offset % 8)
}

fn get_wval_limb_mine(scalar: &Scalar256, offset: usize, bits: usize) -> LimbT {
    let bytes = scalar.as_u8();
    let mut head = offset / 8;

    let top = (offset + bits - 1) / 8 - offset / 8 + 1;

    let mut out = bytes[head] as LimbT;
    for i in 1..4 {
        if i < top {
            head += 1;
            out |= (bytes[head] as LimbT) << (i * 8);
        }
    }

    out >> (offset % 8)
}

fn accumulate_wrapper<TG1: G1 + G1GetFp<TG1Fp>, TG1Fp: G1Fp, TG1Affine: G1Affine<TG1, TG1Fp>>(
    sum: *mut TG1,
    points: *mut TG1,
    n: usize,
) {
    unsafe {
        accumulate::<TG1, TG1Fp, TG1Affine>(sum, points, n);
    }
}

unsafe fn accumulate<TG1: G1 + G1GetFp<TG1Fp>, TG1Fp: G1Fp, TG1Affine: G1Affine<TG1, TG1Fp>>(
    sum: *mut TG1,
    mut points: *mut TG1,
    mut n: usize,
) {
    // ptype *dst; \
    let mut dst;
    // void *mul_acc; \
    let mut mul_acc;
    // size_t i; \
    let mut i;

    // while (n >= 16) { \
    while n >= 16 {
        // if (n & 1) \
        if n & 1 == 1 {
            // ptype##_dadd_affine(sum, sum, (const ptype##_affine *)points++); \
            dadd_affine_ptrs(sum, sum, points as *const TG1Affine);
            points = points.add(1);
        }
        // n /= 2; \
        n /= 2;
        // for (mul_acc = NULL, i = n; i--; mul_acc = points->Z, points += 2) \
        mul_acc = null();
        i = n;
        while i > 0 {
            i -= 1;

            // ptype##_head(points, mul_acc); \
            head::<TG1, TG1Fp, TG1Affine>(points, mul_acc);

            mul_acc = (*points).z() as *const TG1Fp;
            points = points.add(2);
        }

        // reciprocal_##field(points[-2].Z, points[-2].Z); /* 1/∏ Zi */ \
        *(*points.sub(2)).z_mut() = (*points.sub(2)).z().inverse().unwrap();

        // for (dst = points, i = n; --i;) { \
        dst = points;
        i = n;
        while i > 1 {
            i -= 1;
            // dst--; points -= 2; \
            dst = dst.sub(1);
            points = points.sub(2);
            // mul_##field(points[-2].Z, points[0].Z, points[-2].Z); \
            *(*points.sub(2)).z_mut() = (*points).z().mul_fp((*points.sub(2)).z());
            // ptype##_tail(dst, points, points[-2].Z); \
            tail::<TG1, TG1Fp, TG1Affine>(dst, points, (*points.sub(2)).z_mut() as *mut TG1Fp);
            // mul_##field(points[-2].Z, points[0].Z, points[1].Z); \
            *(*points.sub(2)).z_mut() = (*points).z().mul_fp((*points.add(1)).z());
        }
        // dst--; points -= 2; \
        dst = dst.sub(1);
        points = points.sub(2);
        // ptype##_tail(dst, points, points[0].Z); \
        tail::<TG1, TG1Fp, TG1Affine>(dst, points, (*points).z_mut() as *mut TG1Fp);
        // points = dst; \
        points = dst;
    }
    // while (n--) \
    while n > 0 {
        n -= 1;
        dadd_affine_ptrs(sum, sum, points as *const TG1Affine);
        points = points.add(1);
        //     ptype##_dadd_affine(sum, sum, (const ptype##_affine *)points++); \
    }
}

struct AddOrDbl<TG1Fp: G1Fp> {
    h: TG1Fp,
    r: TG1Fp,
    sx: TG1Fp,
}

unsafe fn dadd_affine_ptrs<
    TG1: G1 + G1GetFp<TG1Fp>,
    TG1Fp: G1Fp,
    TG1Affine: G1Affine<TG1, TG1Fp>,
>(
    out: *mut TG1,
    p1p: *const TG1,
    p2p: *const TG1Affine,
) {
    // if true {
    //     *out = dadd_affine(&*p1p, &*p2p);
    //     return;
    // }

    let mut p3x;
    let mut p3y;
    let mut p3z;
    let mut add = AddOrDbl {
        h: TG1Fp::default(),
        r: TG1Fp::default(),
        sx: TG1Fp::default(),
    };
    let mut dbl = AddOrDbl {
        h: TG1Fp::default(),
        r: TG1Fp::default(),
        sx: TG1Fp::default(),
    };

    let p1inf;
    let p2inf;
    let is_dbl;

    let p1 = &*p1p;
    let p2 = &*p2p;

    p2inf = p2.is_zero();
    dbl.sx = p2.x().add_fp(p2.x());
    dbl.r = p2.x().square();
    dbl.r = {
        let tmp = dbl.r.add_fp(&dbl.r);
        dbl.r.add_fp(&tmp)
    };
    dbl.h = p2.y().add_fp(p2.y());

    p1inf = p1.z().is_zero();
    add.h = p1.z().square(); /* Z1^2 */
    add.r = add.h.mul_fp(p1.z()); /* Z1^3 */
    add.r = add.r.mul_fp(p2.y()); /* S2 = Y2*Z1^3 */
    add.r = add.r.sub_fp(p1.y()); /* R = S2-Y1 */

    add.h = add.h.mul_fp(p2.x()); /* U2 = X2*Z1^2 */

    add.sx = add.h.add_fp(p1.x()); /* sx = X1+U2 */
    add.h = add.h.sub_fp(p1.x()); /* H = U2-X1 */

    p3z = add.h.mul_fp(p1.z()); /* Z3 = H*Z1 */

    // /* make the choice between addition and doubling */
    is_dbl = add.h.is_zero() && add.r.is_zero();

    if is_dbl {
        p3x = *p2.x();
        p3y = *p2.y();
        p3z = dbl.h;
        add.h = dbl.h;
        add.r = dbl.r;
        add.sx = dbl.sx;
    } else {
        p3x = *p1.x();
        p3y = *p1.y();
    }
    /* |p3| and |add| hold all inputs now, |p3| will hold output */

    dbl.h = add.h.square(); /* H^2 */
    dbl.r = dbl.h.mul_fp(&add.h); /* H^3 */
    dbl.r = dbl.r.mul_fp(&p3y); /* H^3*S1 */
    p3y = dbl.h.mul_fp(&p3x); /* H^2*U1 */

    dbl.h = dbl.h.mul_fp(&add.sx); /* H^2*sx */
    p3x = add.r.square(); /* R^2 */
    p3x = p3x.sub_fp(&dbl.h); /* X3 = R^2-H^2*sx */

    p3y = p3y.sub_fp(&p3x); /* H^2*U1-X3 */
    p3y = p3y.mul_fp(&add.r); /* R*(H^2*U1-X3) */
    p3y = p3y.sub_fp(&dbl.r); /* Y3 = R*(H^2*U1-X3)-H^3*S1 */

    if p1inf {
        p3x = *p2.x();
        p3y = *p2.y();
        p3z = TG1Fp::one();
    }
    let outr = &mut *out;
    if p2inf {
        *outr.x_mut() = *p1.x();
        *outr.y_mut() = *p1.y();
        *outr.z_mut() = *p1.z();
    } else {
        *outr.x_mut() = p3x;
        *outr.y_mut() = p3y;
        *outr.z_mut() = p3z;
    }
}

unsafe fn head<TG1: G1 + G1GetFp<TG1Fp>, TG1Fp: G1Fp, TG1Affine: G1Affine<TG1, TG1Fp>>(
    ab: *mut TG1,
    mul_acc: *const TG1Fp,
) {
    // ptype *A = AB, *B = AB+1; \
    let mut a = &mut *ab;
    let mut b = &mut *ab.add(1);
    // limb_t inf = vec_is_zero(A, sizeof(ptype##_affine)) | \
    //              vec_is_zero(B, sizeof(ptype##_affine));  \
    let mut inf = (a.x().is_zero() && a.y().is_zero()) || (b.x().is_zero() && b.y().is_zero());
    // static const vec##bits zero = { 0 }; \
    let zero = TG1Fp::zero();

    // sub_##field(B->Z, B->X, A->X);		/* X2-X1  */ \
    *b.z_mut() = b.x().sub_fp(a.x());
    // add_##field(B->X, B->X, A->X);		/* X2+X1  */ \
    *b.x_mut() = b.x().add_fp(a.x());
    // add_##field(A->Z, B->Y, A->Y);		/* Y2+Y1  */ \
    *a.z_mut() = b.y().add_fp(a.y());
    // sub_##field(B->Y, B->Y, A->Y);		/* Y2-Y1  */ \
    *b.y_mut() = b.y().sub_fp(a.y());
    // if (vec_is_zero(B->Z, sizeof(B->Z))) {	/* X2==X1 */ \
    if b.z().is_zero() {
        // inf = vec_is_zero(A->Z, sizeof(A->Z));	\
        inf = a.z().is_zero();
        // vec_select(B->X, A->Z, B->X, sizeof(B->X), inf); \
        if inf {
            *b.x_mut() = *a.z();
        }
        // sqr_##field(B->Y, A->X);		\
        *b.y_mut() = a.x().square();
        // mul_by_3_##field(B->Y, B->Y);		/* 3*X1^2 */ \
        *b.y_mut() = {
            let tmp = b.y().add_fp(b.y());
            tmp.add_fp(b.y())
        };
        // vec_copy(B->Z, A->Z, sizeof(B->Z));	/* 2*Y1   */ \
        *b.z_mut() = *a.z();
    } /* B->Y is numenator    */
    /* B->Z is denominator  */
    // vec_select(A->X, B->X, A->X, sizeof(A->X), inf); \
    // vec_select(A->Y, A->Z, A->Y, sizeof(A->Y), inf); \
    // vec_select(A->Z, one,  B->Z, sizeof(A->Z), inf); \
    // vec_select(B->Z, zero, B->Z, sizeof(B->Z), inf); \
    if inf {
        *a.x_mut() = *b.x();
        *a.y_mut() = *a.z();
        *a.z_mut() = TG1Fp::one();
        *b.z_mut() = TG1Fp::zero();
    } else {
        *a.z_mut() = *b.z();
    }
    // if (mul_acc != NULL) {\
    if (!mul_acc.is_null()) {
        // mul_##field(A->Z, A->Z, *mul_acc);	/* chain multiplication */\
        *a.z_mut() = a.z().mul_fp(&*mul_acc);
    }
}

unsafe fn tail<TG1: G1 + G1GetFp<TG1Fp>, TG1Fp: G1Fp, TG1Affine: G1Affine<TG1, TG1Fp>>(
    d: *mut TG1,
    ab: *mut TG1,
    lambda: *mut TG1Fp,
) {
    // ptype *A = AB, *B = AB+1; \
    let mut a = &mut *ab;
    let mut b = &mut *ab.add(1);
    let mut dr = &mut *d;
    // vec##bits llambda; \
    let mut llambda;
    // limb_t inf = vec_is_zero(B->Z, sizeof(B->Z)); \
    let inf = b.z().is_zero();

    let mut lambdar = &mut *lambda;
    // mul_##field(*lambda, *lambda, B->Y);		/* λ = (Y2-Y1)/(X2-X1)  */ \
    *lambdar = lambdar.mul_fp(b.y());
    /* alt. 3*X1^2/2*Y1     */

    // sqr_##field(llambda, *lambda); \
    llambda = lambdar.square();
    // sub_##field(D->X, llambda, B->X);		/* X3 = λ^2-X1-X2       */ \
    *dr.x_mut() = llambda.sub_fp(b.x());

    // sub_##field(D->Y, A->X, D->X);   \
    *dr.y_mut() = a.x().sub_fp(dr.x());
    // mul_##field(D->Y, D->Y, *lambda); \
    *dr.y_mut() = dr.y().mul_fp(&lambdar);
    // sub_##field(D->Y, D->Y, A->Y);		/* Y3 = λ*(X1-X3)-Y1    */ \
    *dr.y_mut() = dr.y().sub_fp(a.y());

    // vec_select(D->X, A->X, D->X, 2*sizeof(D->X), inf); \
    // vec_select(B->Z, one, B->Z, sizeof(B->Z), inf); \
    if inf {
        *dr.x_mut() = *a.x();
        *dr.y_mut() = *a.y();
        *b.z_mut() = TG1Fp::one();
    }
}

fn booth_encode(wval: LimbT, sz: usize) -> LimbT {
    let mask = 0u64.wrapping_sub(wval >> sz);
    let wval = (wval + 1) >> 1;
    (wval ^ mask).wrapping_sub(mask)
}
