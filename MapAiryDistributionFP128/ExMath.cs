using MultiPrecision;

namespace MapAiryDistributionFP128 {
    internal static class ExMath {

        public static MultiPrecision<Pow2.N4> Square(MultiPrecision<Pow2.N4> x) => x * x;
        public static MultiPrecision<Pow2.N4> Cube(MultiPrecision<Pow2.N4> x) => x * x * x;

        public static MultiPrecision<Pow2.N4> Pow3d2(MultiPrecision<Pow2.N4> x) {
            return MultiPrecision<Pow2.N4>.Sqrt(Cube(x));
        }

        public static MultiPrecision<Pow2.N4> Pow2d3(MultiPrecision<Pow2.N4> x) {
            return MultiPrecision<Pow2.N4>.Cbrt(Square(x));
        }
    }
}