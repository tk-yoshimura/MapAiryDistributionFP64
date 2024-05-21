// Author and Approximation Formula Coefficient Generator: T.Yoshimura
// Github: https://github.com/tk-yoshimura
// Original Code: https://github.com/tk-yoshimura/MapAiryDistributionFP64

using MapAiryDistributionFP64.InternalUtils;
using MapAiryDistributionFP64.RandomGeneration;
using System.Collections.ObjectModel;
using System.Diagnostics;
using static System.Double;

namespace MapAiryDistributionFP64 {
    public class MapAiryDistribution {

        public double Mu { get; }

        public double C { get; }

        private readonly double c_inv;

        private static readonly double mode_base = -1.16158727113597068525;
        private static readonly double median_base = -0.71671068545502205332;
        private static readonly double entropy_base = 2.00727681841065634600;

        public MapAiryDistribution() : this(mu: 0d, c: 1d) { }

        public MapAiryDistribution(double c) : this(mu: 0d, c: c) { }

        public MapAiryDistribution(double mu, double c) {
            if (!IsFinite(mu)) {
                throw new ArgumentOutOfRangeException(nameof(mu), "Invalid location parameter.");
            }
            if (!(c > 0 && IsFinite(c))) {
                throw new ArgumentOutOfRangeException(nameof(c), "Invalid scale parameter.");
            }

            Mu = mu;
            C = c;

            c_inv = 1d / c;
        }

        public double PDF(double x) {
            double u = (x - Mu) * c_inv;

            if (IsNaN(u)) {
                return NaN;
            }
            if (IsInfinity(u)) {
                return 0d;
            }

            double pdf = PDFPade.Value(u) * c_inv;

            return pdf;
        }

        public double CDF(double x, Interval interval = Interval.Lower) {
            double u = (x - Mu) * c_inv;

            if (IsNaN(u)) {
                return NaN;
            }

            double cdf = CDFPade.Value(u, interval != Interval.Lower);

            return cdf;
        }

        public double Quantile(double p, Interval interval = Interval.Lower) {
            if (!(p >= 0d && p <= 1d)) {
                return NaN;
            }

            double x = Mu + C * QuantilePade.Value(p, interval != Interval.Lower);

            return x;
        }

        public double Sample(Random random) {
            double u = random.NextUniformOpenInterval01() - 0.5d;
            double w = random.NextUniformOpenInterval0();

            double cu = CosPi(u);

            double r = -SinPi(u * 1.5d - 0.25d) * Cbrt(2d * Log(w) / (CosPi(u * 0.5d - 0.25d) * cu * cu));
            double v = r * C + Mu;

            return v;
        }

        public double Median => Mu + median_base * C;

        public double Mode => Mu + mode_base * C;

        public double Mean => Mu;

        public double Variance => NaN;

        public double Skewness => NaN;

        public double Kurtosis => NaN;

        public double Entropy => entropy_base + Log(C);

        public double Alpha => 1.5d;

        public double Beta => 1d;

        public static MapAiryDistribution operator +(MapAiryDistribution dist1, MapAiryDistribution dist2) {
            return new(dist1.Mu + dist2.Mu, ExMath.Pow2d3(ExMath.Pow3d2(dist1.C) + ExMath.Pow3d2(dist2.C)));
        }

        public static MapAiryDistribution operator -(MapAiryDistribution dist1, MapAiryDistribution dist2) {
            return new(dist1.Mu - dist2.Mu, ExMath.Pow2d3(ExMath.Pow3d2(dist1.C) + ExMath.Pow3d2(dist2.C)));
        }

        public static MapAiryDistribution operator +(MapAiryDistribution dist, double s) {
            return new(dist.Mu + s, dist.C);
        }

        public static MapAiryDistribution operator -(MapAiryDistribution dist, double s) {
            return new(dist.Mu - s, dist.C);
        }

        public static MapAiryDistribution operator *(MapAiryDistribution dist, double k) {
            return new(dist.Mu * k, dist.C * k);
        }

        public static MapAiryDistribution operator /(MapAiryDistribution dist, double k) {
            return new(dist.Mu / k, dist.C / k);
        }

        public override string ToString() {
            return $"{typeof(MapAiryDistribution).Name}[mu={Mu},c={C}]";
        }

        public string Formula => "p(x; mu, c) := stable_distribution(x; alpha = 3/2, beta = 1, mu, c)";

        private static class PDFPade {
            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_0_1 = new(
                new ReadOnlyCollection<double>([
                    1.97516171847191855610e-1,
                    3.67488253628465083737e-2,
                    -9.73242224038828612673e-4,
                    2.32207514136635673061e-3,
                    5.69067907423210669037e-5,
                    -6.02637387141524535193e-5,
                    1.04960324426666933327e-5,
                    -6.58470237954242016920e-7,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    7.09464351647314165710e-1,
                    3.66413036246461392316e-1,
                    1.10947882302862241488e-1,
                    2.65928486676817177159e-2,
                    3.75507284977386290874e-3,
                    4.03789594641339005785e-4,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_1_2 = new(
                new ReadOnlyCollection<double>([
                    1.06251243013238748252e-1,
                    1.38178831205785069108e-2,
                    4.19280374368049006206e-3,
                    8.54607219684690930289e-4,
                    -7.46881084120928210702e-5,
                    1.47110856483345063335e-5,
                    -1.30090180307471994500e-6,
                    5.24801123304330014713e-8,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    8.10853683888611687140e-1,
                    3.89361261627717143905e-1,
                    1.15124062681082170577e-1,
                    2.38803416611949902468e-2,
                    3.08616898814509065071e-3,
                    2.43760043942846261876e-4,
                    1.34538901435238836768e-6,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_2_4 = new(
                new ReadOnlyCollection<double>([
                    5.33842514891989443409e-2,
                    1.23301980674903270971e-2,
                    3.45717831433988631923e-3,
                    3.27034449923176875761e-4,
                    1.20406794831890291348e-5,
                    5.77489170397965604669e-7,
                    -1.15255267205685159063e-7,
                    9.15896323073109992939e-9,
                    -3.14068002815368247985e-10,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    9.08772985520393226044e-1,
                    4.26418573702560818267e-1,
                    1.22033746594868893316e-1,
                    2.27934009200310243172e-2,
                    2.60658999011198623962e-3,
                    1.54461660261435227768e-4,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_4_8 = new(
                new ReadOnlyCollection<double>([
                    1.58950538583133457384e-2,
                    7.47835440063141601948e-3,
                    1.81137244353261478410e-3,
                    2.26935565382135588558e-4,
                    1.43877113825683795505e-5,
                    2.08242747557417233626e-7,
                    -1.54976465724771282989e-9,
                    1.30762989300333026019e-11,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    9.95505437381674174441e-1,
                    4.58882737262511297099e-1,
                    1.25031310192148865496e-1,
                    2.15727229249904102247e-2,
                    2.33597081566665672569e-3,
                    1.45198998318300328562e-4,
                    3.87962234445835345676e-6,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_8_16 = new(
                new ReadOnlyCollection<double>([
                    3.22517551525042172428e-3,
                    1.12822806030796339659e-3,
                    1.54489389961322571031e-4,
                    9.28479992527909796427e-6,
                    2.06168350199745832262e-7,
                    9.05110751997021418539e-10,
                    -2.15498112371756202097e-12,
                    6.41838355699777435924e-15,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    6.53390465399680164234e-1,
                    1.82759048270449018482e-1,
                    2.80407546367978533849e-2,
                    2.50853443923476718145e-3,
                    1.27671852825846245421e-4,
                    3.28380135691060279203e-6,
                    3.06545317089055335742e-8,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_16_32 = new(
                new ReadOnlyCollection<double>([
                    5.82527663232857270992e-4,
                    6.89502117025124630567e-5,
                    2.24909795087265741433e-6,
                    2.18576787334972903790e-8,
                    3.39014723444178274435e-11,
                    -9.74481309265612390297e-15,
                    -1.13308546492906818388e-16,
                    5.32472028720777735712e-19,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    2.74018883667663396766e-1,
                    2.95901195665990089660e-2,
                    1.57901733512147920251e-3,
                    4.24965124147621236633e-5,
                    5.17522027193205842016e-7,
                    2.00522219276570039934e-9,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_32_64 = new(
                new ReadOnlyCollection<double>([
                    1.03264853379349880039e-4,
                    5.35256306644392405447e-6,
                    9.00657716972118816692e-8,
                    5.34913574042209793720e-10,
                    6.70752605041678779380e-13,
                    -5.30089923101856817552e-16,
                    7.28133811621687143754e-19,
                    -7.38047553655951666420e-22,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    1.29920843258164337377e-1,
                    6.75018577147646502386e-3,
                    1.77694968039695671819e-4,
                    2.46428299911920942946e-6,
                    1.67165053157990942546e-8,
                    4.19496974141131087116e-11,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_limit = new(
                new ReadOnlyCollection<double>([
                    5.98413420602149016910e-1,
                    3.14584075817417883086e-5,
                    1.62977928311793051895e1,
                    -4.12903117172994371875e-4,
                    -1.06404478702135751872e2,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    5.25696892802060720079e-5,
                    4.03600055498020483920e1,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_minus_1_0 = new(
                new ReadOnlyCollection<double>([
                    2.76859868856746781256e-1,
                    1.10489814676299003241e-1,
                    -6.25690643488236678667e-3,
                    -1.17905420222527577236e-3,
                    1.27188963720084274122e-3,
                    -7.20575105181207907889e-5,
                    -2.22575633858411851032e-5,
                    2.94270091008508492304e-6,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    4.98673671503410894284e-1,
                    3.15907666864554716291e-1,
                    8.34463558393629855977e-2,
                    2.71804643993972494173e-2,
                    3.52187050938036578406e-3,
                    7.03072974279509263844e-4,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_minus_2_1 = new(
                new ReadOnlyCollection<double>([
                    2.14483832832989822788e-1,
                    3.72789690317712876663e-1,
                    1.86473650057086284496e-1,
                    1.31182724166379598907e-2,
                    -9.00695064809774432392e-3,
                    3.46884420664996747052e-4,
                    4.88651392754189961173e-4,
                    -6.13516242712196835055e-5,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    1.06478618107122200489e0,
                    4.08809060854459518663e-1,
                    2.66617598099501800866e-1,
                    4.53526315786051807494e-2,
                    2.44078693689626940834e-2,
                    1.52822572478697831870e-3,
                    8.69480001029742502197e-4,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_minus_2_4 = new(
                new ReadOnlyCollection<double>([
                    2.74308494787955998605e-1,
                    4.87765991440983416392e-1,
                    3.84524365110270427617e-1,
                    1.77409497505926097339e-1,
                    5.25612864287310961520e-2,
                    1.01528615034079765421e-2,
                    1.20417225696161842090e-3,
                    6.97462693097107007719e-5,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    1.81256903248465876424e0,
                    1.43959302060852067876e0,
                    6.65882284117861804351e-1,
                    1.97537712781845593211e-1,
                    3.81732970028510912201e-2,
                    4.52767489928026542226e-3,
                    2.62240194911920120003e-4,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_minus_4_8 = new(
                new ReadOnlyCollection<double>([
                    2.67391547707456587286e-1,
                    3.39319035621314371924e-1,
                    1.85434799940724207230e-1,
                    5.63667456320679857693e-2,
                    1.01231164548944177474e-2,
                    1.02501575174439362864e-3,
                    4.60769537123286016400e-5,
                    -4.92754650783224582641e-13,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    1.27271216837333318516e0,
                    6.96551952883867277759e-1,
                    2.11871363524516350422e-1,
                    3.80622887806509632537e-2,
                    3.85400280812991562328e-3,
                    1.73246593953823694311e-4,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_minus_8_16 = new(
                new ReadOnlyCollection<double>([
                    2.66153901932100301337e-1,
                    1.65767350677458230714e-1,
                    4.19801402197670061146e-2,
                    5.39337995172784579558e-3,
                    3.50811247702301287586e-4,
                    9.21758454778883157515e-6,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    6.23092941554668369107e-1,
                    1.57829914506366827914e-1,
                    2.02787979758160988615e-2,
                    1.31903008994475216511e-3,
                    3.46575870637847438219e-5,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_minus_16_32 = new(
                new ReadOnlyCollection<double>([
                    2.65985830928929730672e-1,
                    7.12399432979613322705e-2,
                    7.12711058905939981164e-3,
                    3.15786968248685705045e-4,
                    5.22817461604216528366e-6,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    2.67850712372594664473e-1,
                    2.67975253839917164432e-2,
                    1.18734081496856828219e-3,
                    1.96576354766834858479e-5,
                ])
            );

            private static double PlusValue(double x) {
                Debug.Assert(x >= 0);

                double y;
                if (x <= 1d) {
                    Debug.WriteLine("pade minimum segment passed");

                    y = ApproxUtil.Pade(x, pade_plus_0_1);
                }
                else if (x <= 2d) {
                    y = ApproxUtil.Pade(x - 1d, pade_plus_1_2);
                }
                else if (x <= 4d) {
                    y = ApproxUtil.Pade(x - 2d, pade_plus_2_4);
                }
                else if (x <= 8d) {
                    y = ApproxUtil.Pade(x - 4d, pade_plus_4_8);
                }
                else if (x <= 16d) {
                    y = ApproxUtil.Pade(x - 8d, pade_plus_8_16);
                }
                else if (x <= 32d) {
                    y = ApproxUtil.Pade(x - 16d, pade_plus_16_32);
                }
                else if (x <= 64d) {
                    y = ApproxUtil.Pade(x - 32d, pade_plus_32_64);
                }
                else {
                    double u = 1d / ExMath.Pow3d2(x);

                    y = ApproxUtil.Pade(u, pade_plus_limit) * u / x;
                }

                return y;
            }

            private static double MinusValue(double x) {
                Debug.Assert(x <= 0);

                x = -x;

                if (x <= 2d) {
                    double y;
                    if (x <= 1d) {
                        y = ApproxUtil.Pade(1d - x, pade_minus_1_0);
                    }
                    else {
                        y = ApproxUtil.Pade(2d - x, pade_minus_2_1);
                    }

                    return y;
                }
                else if (x <= 32d) {
                    double v;
                    if (x <= 4d) {
                        v = ApproxUtil.Pade(x - 2d, pade_minus_2_4);
                    }
                    else if (x <= 8d) {
                        v = ApproxUtil.Pade(x - 4d, pade_minus_4_8);
                    }
                    else if (x <= 16d) {
                        v = ApproxUtil.Pade(x - 8d, pade_minus_8_16);
                    }
                    else {
                        v = ApproxUtil.Pade(x - 16d, pade_minus_16_32);
                    }

                    double y = v * Sqrt(x) * Exp(-(2d * x * x * x) / 27d);

                    return y;
                }
                else {
                    return 0d;
                }
            }

            public static double Value(double x) {
                return (x >= 0d) ? PlusValue(x) : (x <= 0d) ? MinusValue(x) : NaN;
            }
        }

        private static class CDFPade {
            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_0_1 = new(
                new ReadOnlyCollection<double>([
                    3.33333333333333333333e-1,
                    7.49532137610545010591e-2,
                    9.25326921848155048716e-3,
                    6.59133092365796208900e-3,
                    -5.21942678326323374113e-4,
                    8.22766804917461941348e-5,
                    -3.97941251650023182117e-6,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    8.17408156824742736411e-1,
                    3.57041011418415988268e-1,
                    1.04580353775369716002e-1,
                    1.87521616934129432292e-2,
                    2.33232161135637085535e-3,
                    7.31285352607895467310e-5,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_1_2 = new(
                new ReadOnlyCollection<double>([
                    1.84196970581015939888e-1,
                    -1.19398028299089933853e-3,
                    1.21954054797949597854e-2,
                    -9.37912675685073154845e-4,
                    1.66651954077980453212e-4,
                    -1.33271812303025233648e-5,
                    5.35982226125013888796e-7,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    5.70352826101668448273e-1,
                    1.98852010141232271304e-1,
                    3.64864882318453496161e-2,
                    4.22173125405065522298e-3,
                    1.20079284386796600356e-4,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_2_4 = new(
                new ReadOnlyCollection<double>([
                    1.07409273397524124098e-1,
                    3.83900318969331880402e-2,
                    1.17926652359826576790e-2,
                    1.52181625871479030046e-3,
                    1.50703424417132565662e-4,
                    2.10117959279448106308e-6,
                    1.97360985832285866640e-8,
                    -1.06076300080048408251e-9,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    8.54435380513870673497e-1,
                    3.66021233157880878411e-1,
                    9.42985570806905160687e-2,
                    1.54122343653998564507e-2,
                    1.49849056258932455548e-3,
                    6.94290406268856211707e-5,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_4_8 = new(
                new ReadOnlyCollection<double>([
                    4.70720199535228802538e-2,
                    2.67200763833749070079e-2,
                    7.37400551855064729769e-3,
                    1.10592441765001623699e-3,
                    9.15846028547400212588e-5,
                    3.17801522553862136789e-6,
                    2.03102753319827713542e-8,
                    -5.16172854149066643529e-11,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    9.05317644829451086870e-1,
                    3.73713496637025562492e-1,
                    8.94434672792094976627e-2,
                    1.31846542255347106087e-2,
                    1.16680596342421447100e-3,
                    5.44719256441278863300e-5,
                    8.73131209154185067287e-7,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_8_16 = new(
                new ReadOnlyCollection<double>([
                    1.74847564444513000450e-2,
                    6.00209162595027323742e-3,
                    7.86550260761375576075e-4,
                    4.46682547335758521734e-5,
                    9.51329761417139273391e-7,
                    4.10313065114362712333e-9,
                    -9.81286503831545640189e-12,
                    2.98763969872672156104e-14,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    5.27732094554221674504e-1,
                    1.14330643482604301178e-1,
                    1.27722341942374066265e-2,
                    7.54563340152441778517e-4,
                    2.13377039814057925832e-5,
                    2.09670987094350618690e-7,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_16_32 = new(
                new ReadOnlyCollection<double>([
                    6.22684103170563193015e-3,
                    1.34714356588780958096e-3,
                    9.51289465377874891896e-5,
                    2.64918464474843134081e-6,
                    2.66703857491046795285e-8,
                    5.42037888457985833156e-11,
                    -6.18017115447736427379e-14,
                    9.11626234402148561268e-17,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    3.09895694991285975774e-1,
                    3.69874670435930773471e-2,
                    2.15708854325146400153e-3,
                    6.35345408451056881884e-5,
                    8.65722805575670770555e-7,
                    4.03153189557220023202e-9,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_32_64 = new(
                new ReadOnlyCollection<double>([
                    2.20357145727036120652e-3,
                    1.45412555771401325111e-4,
                    3.27819006009093198652e-6,
                    2.96786786716623870006e-8,
                    9.54192199129339742308e-11,
                    5.71421706870777687254e-14,
                    -1.48321866072033823195e-17,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    1.12851983233980279746e-1,
                    4.94650928817638043712e-3,
                    1.05447405092956497114e-4,
                    1.11578464291338271178e-6,
                    5.27522295397347842625e-9,
                    7.95786524903707645399e-12,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_limit = new(
                new ReadOnlyCollection<double>([
                    3.98942280401432677940e-1,
                    2.89752186412133782995e-2,
                    4.67360459917040710474e0,
                    -1.26770824563800250704e-1,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    7.26301023103568827709e-2,
                    1.60899894281099149848e1,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_minus_1_0 = new(
                new ReadOnlyCollection<double>([
                    4.23238998449671083670e-1,
                    4.95353582976475183891e-1,
                    2.45823281826037784270e-1,
                    7.29726507468813920788e-2,
                    1.63332856186819713346e-2,
                    2.82514634871307516142e-3,
                    2.66220579589280704089e-4,
                    3.09442180091323751049e-6,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    5.16241922223786900600e-1,
                    2.75690727171711638879e-1,
                    7.18707184893542884080e-2,
                    1.87136800286819336797e-2,
                    2.38383441176345054929e-3,
                    3.23509126477812051983e-4,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_minus_2_1 = new(
                new ReadOnlyCollection<double>([
                    1.62598955251978523175e-1,
                    2.30154661502402196205e-1,
                    1.29233975368291684522e-1,
                    3.80919553916980965587e-2,
                    8.17724414618808505948e-3,
                    1.95816800210481122544e-3,
                    3.35259917978421935141e-4,
                    1.22071311320012805777e-5,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    9.63771793313770952352e-2,
                    2.23602260938227310054e-1,
                    9.21944797677283179038e-3,
                    1.82181136341939651516e-2,
                    1.11216849284965970458e-4,
                    5.57446347676836375810e-4,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_minus_2_4 = new(
                new ReadOnlyCollection<double>([
                    5.88176189476056502705e-1,
                    6.02287088109671443912e-1,
                    2.57625533207709555766e-1,
                    7.18618327959270311402e-2,
                    1.25253508822578071586e-2,
                    1.38226766788726569279e-3,
                    7.32620930807726458318e-5,
                    6.03054572349954521563e-7,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    9.54199287706410202241e-1,
                    4.75472282467470181184e-1,
                    1.46535255063482444037e-1,
                    3.01738952830506685803e-2,
                    3.99499614754061850707e-3,
                    3.12152806855533452870e-4,
                    8.05172630869950926085e-6,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_minus_4_8 = new(
                new ReadOnlyCollection<double>([
                    5.51472757643673529728e-1,
                    4.31610102565728326076e-1,
                    1.48917299048993249813e-1,
                    2.91436404351984593223e-2,
                    3.28882396835475166496e-3,
                    1.98804328328017984553e-4,
                    4.64616080021023007510e-6,
                    1.56438181721678316854e-8,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    8.56039144443074807236e-1,
                    3.34960640112789039792e-1,
                    7.52361661319904780994e-2,
                    1.02014114984162352488e-2,
                    7.97578941836281798696e-4,
                    3.03148087501011803869e-5,
                    3.26204774056531450781e-7,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_minus_8_16 = new(
                new ReadOnlyCollection<double>([
                    4.18065737742332603636e-1,
                    2.01623556896170223584e-1,
                    3.95777759442922635350e-2,
                    3.93741411702090093612e-3,
                    2.00425740936518942152e-4,
                    4.65171940494245225320e-6,
                    3.71059603439031417668e-8,
                    4.25509635435515655505e-11,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    5.40385920845795746296e-1,
                    1.21525463547886545439e-1,
                    1.43424513283784662739e-2,
                    9.21390853863369654347e-4,
                    3.01133286018215339596e-5,
                    4.14450918250118242069e-7,
                    1.51212616739642117089e-9,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_minus_16_32 = new(
                new ReadOnlyCollection<double>([
                    2.98743266203527920889e-1,
                    4.83180073178398517903e-2,
                    2.69520487031514045597e-3,
                    5.91803613302232034530e-5,
                    4.28395710187165863411e-7,
                    4.53343398774819924776e-10,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    1.92698084903936492789e-1,
                    1.35679707446151988505e-2,
                    4.15604984260210557217e-4,
                    5.08716330609178693652e-6,
                    1.66933196369632373584e-8,
                ])
            );

            private static double PlusValue(double x) {
                Debug.Assert(x >= 0);

                double y;
                if (x <= 1d) {
                    Debug.WriteLine("pade minimum segment passed");

                    y = ApproxUtil.Pade(x, pade_plus_0_1);
                }
                else if (x <= 2d) {
                    y = ApproxUtil.Pade(x - 1d, pade_plus_1_2);
                }
                else if (x <= 4d) {
                    y = ApproxUtil.Pade(x - 2d, pade_plus_2_4);
                }
                else if (x <= 8d) {
                    y = ApproxUtil.Pade(x - 4d, pade_plus_4_8);
                }
                else if (x <= 16d) {
                    y = ApproxUtil.Pade(x - 8d, pade_plus_8_16);
                }
                else if (x <= 32d) {
                    y = ApproxUtil.Pade(x - 16d, pade_plus_16_32);
                }
                else if (x <= 64d) {
                    y = ApproxUtil.Pade(x - 32d, pade_plus_32_64);
                }
                else {
                    double u = 1d / ExMath.Pow3d2(x);

                    y = ApproxUtil.Pade(u, pade_plus_limit) * u;
                }

                return y;
            }

            private static double MinusValue(double x) {
                Debug.Assert(x <= 0);

                x = -x;

                if (x <= 2d) {
                    double y;
                    if (x <= 1d) {
                        Debug.WriteLine("pade minimum segment passed");

                        y = ApproxUtil.Pade(1d - x, pade_minus_1_0);
                    }
                    else {
                        y = ApproxUtil.Pade(2d - x, pade_minus_2_1);
                    }

                    return y;
                }
                else if (x <= 32d) {
                    double v;
                    if (x <= 4d) {
                        v = ApproxUtil.Pade(x - 2d, pade_minus_2_4);
                    }
                    else if (x <= 8d) {
                        v = ApproxUtil.Pade(x - 4d, pade_minus_4_8);
                    }
                    else if (x <= 16d) {
                        v = ApproxUtil.Pade(x - 8d, pade_minus_8_16);
                    }
                    else {
                        v = ApproxUtil.Pade(x - 16d, pade_minus_16_32);
                    }

                    double y = v * Exp(-(2d * x * x * x) / 27d) / x;

                    return y;
                }
                else {
                    return 0d;
                }
            }

            public static double Value(double x, bool complementary) {
                if (x >= 0d) {
                    return complementary ? PlusValue(x) : 1d - PlusValue(x);
                }
                else if (x <= 0d) {
                    return complementary ? 1d - MinusValue(x) : MinusValue(x);
                }
                else {
                    return NaN;
                }
            }
        }

        private static class QuantilePade {
            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_upper_0p125_0p25 = new(
                new ReadOnlyCollection<double>([
                    1.70276979914029733585e0,
                    2.09991992116646276165e1,
                    2.26775403775298867998e1,
                    -4.85384304722129472833e2,
                    -1.47107146466495573999e3,
                    -7.08748473959943943929e1,
                    1.54245210917147215257e3,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    2.13092357122115486375e1,
                    1.57318281834689144053e2,
                    4.42261730187813035957e2,
                    2.10814431586717588454e2,
                    -6.36700983439599552504e2,
                    -2.82923881266630617596e2,
                    1.36613971025062750340e2,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_upper_0p25_0p5 = new(
                new ReadOnlyCollection<double>([
                    4.81512108276093785320e-1,
                    -2.74296316128959647914e0,
                    -3.29973875964825685757e1,
                    -4.87536980816224603581e1,
                    8.22233203036734027999e1,
                    1.21654607908452130093e2,
                    -6.66681853240657307279e1,
                    -4.28101952511581488588e1,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    8.20189490825315245036e0,
                    1.63469912146101848441e1,
                    -1.52740920318273920072e1,
                    -5.41684560257839409762e1,
                    6.51733677169299416471e0,
                    3.93092001388102589237e1,
                    -9.59983666140749481195e-1,
                    -9.95648827557655863699e-1,
                    -1.32007124426778083829e0,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_upper_expm3_4 = new(
                new ReadOnlyCollection<double>([
                    4.25692449785074345588e-1,
                    3.10963501706596356267e-1,
                    2.91357806215297069863e-2,
                    2.34716342676849303244e-2,
                    5.83137296293361915583e-3,
                    3.71792415497884868748e-4,
                    1.59538372221030642757e-4,
                    4.74040834029330213692e-6,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    4.14801234415100707213e-1,
                    1.04693730144480856638e-1,
                    3.81581484862997435076e-2,
                    8.95334009127358617362e-3,
                    1.43316686981760147226e-3,
                    1.81367766024620080990e-4,
                    1.54779999748286671973e-5,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_upper_expm4_8 = new(
                new ReadOnlyCollection<double>([
                    5.07341098045260541890e-1,
                    3.11771145411143166935e-1,
                    1.74515601081894060888e-1,
                    8.46576990174024231338e-2,
                    2.57510090204322149315e-2,
                    8.26605326867021684811e-3,
                    1.73081423934722046819e-3,
                    3.36314161099011673569e-4,
                    4.50990441180388912803e-5,
                    4.53513191985642134268e-6,
                    2.62304611053075404923e-7,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    5.28225379952156944029e-1,
                    3.49662079845715371907e-1,
                    1.45408903426879603625e-1,
                    5.06773501409016231879e-2,
                    1.45385556714043243731e-2,
                    3.31235831325018043744e-3,
                    6.06977554525543056050e-4,
                    8.42406730405209749492e-5,
                    8.32337989541696717905e-6,
                    4.84923196546857128337e-7,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_upper_expm8_16 = new(
                new ReadOnlyCollection<double>([
                    5.41774626094491510395e-1,
                    4.11060141334529017898e-1,
                    1.48195601801946264526e-1,
                    3.33881552814492855873e-2,
                    5.20893974732203890418e-3,
                    5.84734765774178832854e-4,
                    4.71028150898133935445e-5,
                    2.59185739450631464618e-6,
                    7.77428184258777394627e-8,
                    2.51255632629650930196e-14,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    7.58341767924960527280e-1,
                    2.73511775500642961539e-1,
                    6.16011987856129890130e-2,
                    9.61296002312356116021e-3,
                    1.07890675777726076554e-3,
                    8.69223632953458271977e-5,
                    4.78248875031756169279e-6,
                    1.43460852065144859304e-7,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_upper_expm16_32 = new(
                new ReadOnlyCollection<double>([
                    5.41926067826974905066e-1,
                    4.86926556246548518715e-1,
                    2.11963908288176005856e-1,
                    5.92200639925655576883e-2,
                    1.18859816815542567438e-2,
                    1.76833662992855443754e-3,
                    2.21226152157950219596e-4,
                    1.50444847316426133872e-5,
                    1.87458213915373906356e-6,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    8.98511036742503939380e-1,
                    3.91130673008184655152e-1,
                    1.09277016228474605069e-1,
                    2.19328471889880028208e-2,
                    3.26305879571349016107e-3,
                    4.08222014684743492069e-4,
                    2.77611385768697969181e-5,
                    3.45911046256304795257e-6,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_upper_expm32_48 = new(
                new ReadOnlyCollection<double>([
                    5.41926070139289008291e-1,
                    6.93835278521566240557e-1,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    1.28031352753196139513e0,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_lower_0p125_0p25 = new(
                new ReadOnlyCollection<double>([
                    -2.18765177572396469657e0,
                    -3.65752788934974426531e1,
                    -1.81144810822028903904e2,
                    -1.22434531262312950288e2,
                    8.99451018491165823831e2,
                    9.11333307522308410858e2,
                    -8.76285742384616909177e2,
                    -2.33786726970025938837e2,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    1.91797638291395345792e1,
                    1.24293724082506952768e2,
                    2.82393116012902543276e2,
                    -1.80472369158936285558e1,
                    -5.31764390192922827093e2,
                    -5.60586018315854885788e1,
                    1.21284324755968033098e2,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_lower_0p25_0p375 = new(
                new ReadOnlyCollection<double>([
                    -1.63281240925531302762e0,
                    -4.92351310795930780147e0,
                    1.43448529253101759409e1,
                    3.33182629948094299473e1,
                    -3.06679026539368582747e1,
                    -2.87298447423841965301e1,
                    1.31575930750093554120e1,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    5.38761652244702318296e0,
                    2.40932080746189543284e0,
                    -1.69465870062123632126e1,
                    -6.39998944283654848809e0,
                    1.27168434054332272391e1,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_lower_0p375_0p5 = new(
                new ReadOnlyCollection<double>([
                    -1.17326074020471664075e0,
                    1.51461298154568349598e0,
                    1.19979368094343490487e1,
                    -5.94882121521324108164e0,
                    -2.20619749774447254528e1,
                    7.17766543775229176131e0,
                    4.79284243496552841508e0,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    1.76268072706610602584e0,
                    -4.88492535243404839734e0,
                    -5.67524172432687656881e0,
                    6.83327389947131710596e0,
                    2.91338085774159042709e0,
                    -1.41108918944159283950e0,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_lower_expm3_4 = new(
                new ReadOnlyCollection<double>([
                    -2.18765177572396470773e0,
                    -2.19887766409334094428e0,
                    -7.77080107207360785208e-1,
                    -1.15551765136654549650e-1,
                    -6.64711321022529990367e-3,
                    -9.74212491048543799073e-5,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    7.91919722132624625590e-1,
                    2.17415447268626558639e-1,
                    2.41474762519410575392e-2,
                    9.41084107182696904714e-4,
                    6.65754108797614202364e-6,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_lower_expm4_8 = new(
                new ReadOnlyCollection<double>([
                    -2.59822399410385085335e0,
                    -2.24306757759003016244e0,
                    -7.36208578161752060979e-1,
                    -1.15130762650287391576e-1,
                    -8.77652386123688618995e-3,
                    -2.96358888256575251437e-4,
                    -3.33661282483762192446e-6,
                    -4.19292241201527861927e-9,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    7.23065798041556418844e-1,
                    1.96731305131315877264e-1,
                    2.49952034298034383781e-2,
                    1.49149568322111062242e-3,
                    3.66010398525593921460e-5,
                    2.46857713549279930857e-7,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_lower_expm8_16 = new(
                new ReadOnlyCollection<double>([
                    -3.67354365380697580447e0,
                    -1.52181685844845957618e0,
                    -2.40883948836320845233e-1,
                    -1.82424079258401987512e-2,
                    -6.75844978572417703979e-4,
                    -1.11273358356809152121e-5,
                    -6.12797605223700996671e-8,
                    -3.78061321691170114390e-11,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    3.57770840766081587688e-1,
                    4.81290550545412209056e-2,
                    3.02079969075162071807e-3,
                    8.89589626547135423615e-5,
                    1.07618717290978464257e-6,
                    3.57383804712249921193e-9,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_lower_expm16_32 = new(
                new ReadOnlyCollection<double>([
                    -4.92187819510636697128e0,
                    -9.94924018698264727979e-1,
                    -7.69914962772717316098e-2,
                    -2.85558010159310978248e-3,
                    -5.19022578720207406789e-5,
                    -4.19975546950263453259e-7,
                    -1.13886013623971006760e-9,
                    -3.46758191090170732580e-13,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    1.77270673840643360017e-1,
                    1.18099604045834575786e-2,
                    3.66889581757166584963e-4,
                    5.34484782554469770841e-6,
                    3.19694601727035291809e-8,
                    5.24649233511937214948e-11,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_lower_expm32_64 = new(
                new ReadOnlyCollection<double>([
                    -6.41443550638291133784e0,
                    -6.38369359780748328332e-1,
                    -2.43420704406734621618e-2,
                    -4.45274771094277987075e-4,
                    -3.99529078051262843241e-6,
                    -1.59758677464731620413e-8,
                    -2.14338367751477432622e-11,
                    -3.23343844538964435927e-15,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    8.79845511272943785289e-2,
                    2.90839059356197474893e-3,
                    4.48172838083912540123e-5,
                    3.23770691025690100895e-7,
                    9.60156044379859908674e-10,
                    7.81134095049301988435e-13,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_lower_expm64_128 = new(
                new ReadOnlyCollection<double>([
                    -8.23500806363233610938e0,
                    -4.05652655284908839003e-1,
                    -7.65978833819859622912e-3,
                    -6.94194676058731901672e-5,
                    -3.08771646223818451436e-7,
                    -6.12443207313641110962e-10,
                    -4.07882839359528825925e-13,
                    -3.05720104049292610799e-17,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    4.37395212065018405474e-2,
                    7.18654254114820140590e-4,
                    5.50371158026951899491e-6,
                    1.97583864365011234715e-8,
                    2.91169706068202431036e-11,
                    1.17716830382540977039e-14,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_lower_expm128_256 = new(
                new ReadOnlyCollection<double>([
                    -1.04845570631944023913e1,
                    -2.56502856165700644836e-1,
                    -2.40615394566347412600e-3,
                    -1.08364601171893250764e-5,
                    -2.39603255140022514289e-8,
                    -2.36344017673944676435e-11,
                    -7.83146284114485675414e-15,
                    -2.92218240202835807955e-19,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    2.17740414929742679904e-2,
                    1.78084231709097280884e-4,
                    6.78870668961146609668e-7,
                    1.21313439060489363960e-9,
                    8.89917934953781122884e-13,
                    1.79115540847944524599e-16,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_lower_expm256_512 = new(
                new ReadOnlyCollection<double>([
                    -1.32865827226175698181e1,
                    -1.61802434199627472010e-1,
                    -7.55642602577784211259e-4,
                    -1.69457608092375302291e-6,
                    -1.86612389867293722402e-9,
                    -9.17015770142364635163e-13,
                    -1.51422473889348610974e-16,
                    -2.81661279271583206526e-21,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    1.08518414679241420227e-2,
                    4.42335224797004486239e-5,
                    8.40387821972524402121e-8,
                    7.48486746424527560620e-11,
                    2.73676810622938942041e-14,
                    2.74588200481263214866e-18,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_lower_expm512_1024 = new(
                new ReadOnlyCollection<double>([
                    -1.67937186583822375593e1,
                    -1.01958138247797604098e-1,
                    -2.37409774265951876695e-4,
                    -2.65483321307104128810e-7,
                    -1.45803536947907216594e-10,
                    -3.57375116523338994342e-14,
                    -2.94401318006358820268e-18,
                    -2.73260616170245224789e-23,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    5.41357843707822974161e-3,
                    1.10082540037527566536e-5,
                    1.04338126042963003178e-8,
                    4.63619608458569600346e-12,
                    8.45781310395535984099e-16,
                    4.23432554226506409568e-20,
                ])
            );

            private static double UpperValue(double x) {
                Debug.Assert(x <= 0.5d);

                if (x >= 0.125d) {
                    double y;
                    if (x <= 0.25d) {
                        y = ApproxUtil.Pade(x - 0.125d, pade_upper_0p125_0p25);
                    }
                    else {
                        y = ApproxUtil.Pade(x - 0.25d, pade_upper_0p25_0p5);
                    }

                    return y;
                }
                else {
                    double v;
                    int exponent = ILogB(x);

                    if (exponent >= -4) {
                        v = ApproxUtil.Pade(-Log2(ScaleB(x, 3)), pade_upper_expm3_4);
                    }
                    else if (exponent >= -8) {
                        v = ApproxUtil.Pade(-Log2(ScaleB(x, 4)), pade_upper_expm4_8);
                    }
                    else if (exponent >= -16) {
                        v = ApproxUtil.Pade(-Log2(ScaleB(x, 8)), pade_upper_expm8_16);
                    }
                    else if (exponent >= -32) {
                        v = ApproxUtil.Pade(-Log2(ScaleB(x, 16)), pade_upper_expm16_32);
                    }
                    else if (exponent >= -48) {
                        v = ApproxUtil.Pade(-Log2(ScaleB(x, 32)), pade_upper_expm32_48);
                    }
                    else {
                        v = 1d / Cbrt(ScaleB(Pi, 1));
                    }

                    double y = v / ExMath.Pow2d3(x);

                    return y;
                }
            }

            private static double LowerValue(double x) {
                Debug.Assert(x <= 0.5d);

                if (x >= 0.125d) {
                    double y;
                    if (x <= 0.25d) {
                        y = ApproxUtil.Pade(x - 0.125d, pade_lower_0p125_0p25);
                    }
                    else if (x <= 0.375d) {
                        y = ApproxUtil.Pade(x - 0.25d, pade_lower_0p25_0p375);
                    }
                    else {
                        y = ApproxUtil.Pade(x - 0.375d, pade_lower_0p375_0p5);
                    }

                    return y;
                }
                else {
                    double y;
                    int exponent = ILogB(x);

                    if (exponent >= -4) {
                        y = ApproxUtil.Pade(-Log2(ScaleB(x, 3)), pade_lower_expm3_4);
                    }
                    else if (exponent >= -8) {
                        y = ApproxUtil.Pade(-Log2(ScaleB(x, 4)), pade_lower_expm4_8);
                    }
                    else if (exponent >= -16) {
                        y = ApproxUtil.Pade(-Log2(ScaleB(x, 8)), pade_lower_expm8_16);
                    }
                    else if (exponent >= -32) {
                        y = ApproxUtil.Pade(-Log2(ScaleB(x, 16)), pade_lower_expm16_32);
                    }
                    else if (exponent >= -64) {
                        y = ApproxUtil.Pade(-Log2(ScaleB(x, 32)), pade_lower_expm32_64);
                    }
                    else if (exponent >= -128) {
                        y = ApproxUtil.Pade(-Log2(ScaleB(x, 64)), pade_lower_expm64_128);
                    }
                    else if (exponent >= -256) {
                        y = ApproxUtil.Pade(-Log2(ScaleB(x, 128)), pade_lower_expm128_256);
                    }
                    else if (exponent >= -512) {
                        y = ApproxUtil.Pade(-Log2(ScaleB(x, 256)), pade_lower_expm256_512);
                    }
                    else if (exponent >= -1024) {
                        y = ApproxUtil.Pade(-Log2(ScaleB(x, 512)), pade_lower_expm512_1024);
                    }
                    else {
                        return NegativeInfinity;
                    }

                    return y;
                }
            }

            public static double Value(double x, bool complementary) {
                if (x > 0.5) {
                    return Value(1d - x, !complementary);
                }

                return complementary ? UpperValue(x) : LowerValue(x);
            }
        }
    }
}