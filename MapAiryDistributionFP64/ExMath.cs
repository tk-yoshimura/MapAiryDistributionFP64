using static System.Double;

namespace MapAiryDistributionFP64 {
    internal static class ExMath {

        public static double Square(double x) => x * x;
        public static double Cube(double x) => x * x * x;

        public static double Pow3d2(double x) {
            if (x == 0d || !IsFinite(x) || int.Abs(double.ILogB((double)x)) < 320) {
                return Sqrt(Cube(x));
            }
            else {
                return Cube(Sqrt(x));
            }
        }

        public static double Pow2d3(double x) {
            if (x == 0d || !IsFinite(x) || int.Abs(double.ILogB((double)x)) < 480) {
                return Cbrt(Square(x));
            }
            else {
                return Square(Cbrt(x));
            }
        }
    }
}