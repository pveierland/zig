// Ported from musl, which is licensed under the MIT license:
// https://git.musl-libc.org/cgit/musl/tree/COPYRIGHT
//
// https://git.musl-libc.org/cgit/musl/tree/src/math/acoshf.c
// https://git.musl-libc.org/cgit/musl/tree/src/math/acosh.c

const std = @import("../std.zig");
const math = std.math;
const expect = std.testing.expect;
const expectApproxEqRel = std.testing.expectApproxEqRel;

/// Returns the hyperbolic arc-cosine of x.
///
/// Special cases:
///  - acosh(x)   = nan if x < 1
///  - acosh(nan) = nan
pub fn acosh(x: anytype) @TypeOf(x) {
    const T = @TypeOf(x);
    return switch (T) {
        f32 => acosh32(x),
        f64 => acosh64(x),
        else => @compileError("acosh not implemented for " ++ @typeName(T)),
    };
}

// acosh(x) = log(x + sqrt(x * x - 1))
fn acosh32(x: f32) f32 {
    const u = @as(u32, @bitCast(x));

    // 0 <= x < 2, invalid if x < 1
    if (u < 0x3F800000 + (1 << 23)) {
        return math.log1p(x - 1 + @sqrt((x - 1) * (x - 1) + 2 * (x - 1)));
    }
    // 2 <= x < 0x1p12
    else if (u < 0x3F800000 + (12 << 23)) {
        return @log(2 * x - 1 / (x + @sqrt(x * x - 1)));
    }
    // x < 0 or x = -0, invalid
    else if (u & 0x80000000 != 0) {
        return math.nan(f32);
    }
    // x >= 0x1p12 or nan
    else {
        return @log(x) + 0.693147180559945309417232121458176568;
    }
}

fn acosh64(x: f64) f64 {
    const u = @as(u64, @bitCast(x));
    const e = u >> 52;

    // 0 <= x < 2, invalid if x < 1
    if (e < 0x3FF + 1) {
        return math.log1p(x - 1 + @sqrt((x - 1) * (x - 1) + 2 * (x - 1)));
    }
    // 2 <= x < 0x1p26
    else if (e < 0x3FF + 26) {
        return @log(2 * x - 1 / (x + @sqrt(x * x - 1)));
    }
    // x < 0 or x = -0, invalid
    else if (e & 0x800 != 0) {
        return math.nan(f64);
    }
    // x >= 0x1p26 or nan
    else {
        return @log(x) + 0.693147180559945309417232121458176568;
    }
}

test acosh {
    try expect(acosh(@as(f32, 1.5)) == acosh32(1.5));
    try expect(acosh(@as(f64, 1.5)) == acosh64(1.5));
}

test acosh32 {
    const epsilon = 0.000001;

    try expect(math.approxEqAbs(f32, acosh32(1.5), 0.962424, epsilon));
    try expect(math.approxEqAbs(f32, acosh32(37.45), 4.315976, epsilon));
    try expect(math.approxEqAbs(f32, acosh32(89.123), 5.183133, epsilon));
    try expect(math.approxEqAbs(f32, acosh32(123123.234375), 12.414088, epsilon));

    // Formula boundaries, with expected values from MPFR 4.2.2
    // (mpfr_acosh, 24 bits, MPFR_RNDN).
    const tolerance = 2 * math.floatEps(f32);
    try expectApproxEqRel(0x1.512428p0, acosh32(2), tolerance);
    try expectApproxEqRel(0x1.205966p3, acosh32(0x1p12), tolerance);
}

test acosh64 {
    const epsilon = 0.000001;

    try expect(math.approxEqAbs(f64, acosh64(1.5), 0.962424, epsilon));
    try expect(math.approxEqAbs(f64, acosh64(37.45), 4.315976, epsilon));
    try expect(math.approxEqAbs(f64, acosh64(89.123), 5.183133, epsilon));
    try expect(math.approxEqAbs(f64, acosh64(123123.234375), 12.414088, epsilon));

    // Formula boundaries, with expected values from MPFR 4.2.2
    // (mpfr_acosh, 53 bits, MPFR_RNDN).
    const tolerance = 2 * math.floatEps(f64);
    try expectApproxEqRel(0x1.5124271980435p0, acosh64(2), tolerance);
    try expectApproxEqRel(0x1.2b708872320e2p4, acosh64(0x1p26), tolerance);
}

test "acosh32.special" {
    try expect(math.isNan(acosh32(math.nan(f32))));
    try expect(math.isNan(acosh32(0.5)));

    for ([_]f32{ 0, -0.0, -1, -1024, -0x1p12, -math.inf(f32) }) |x| {
        try expect(math.isNan(acosh32(x)));
    }

    try expect(math.isPositiveZero(acosh32(1)));
    try expect(math.isPositiveInf(acosh32(math.inf(f32))));
}

test "acosh64.special" {
    try expect(math.isNan(acosh64(math.nan(f64))));
    try expect(math.isNan(acosh64(0.5)));

    // -0x1.8p15 is the input from musl commit 3e80328d.
    for ([_]f64{ 0, -0.0, -1, -0x1.8p15, -0x1p26, -math.inf(f64) }) |x| {
        try expect(math.isNan(acosh64(x)));
    }

    try expect(math.isPositiveZero(acosh64(1)));
    try expect(math.isPositiveInf(acosh64(math.inf(f64))));
}
