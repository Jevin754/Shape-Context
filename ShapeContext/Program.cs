using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ShapeContext
{
    class Program
    {
        static void Main(string[] args)
        {
            Console.WriteLine("H");
        }
    }
    class point
    {
        private double x, y, z;


        public point(double a, double b, double c)
        {
            x = a;
            y = b;
            z = c;
        }

        public double GetX()
        {
            return x;
        }

        public double GetY()
        {
            return y;
        }

        public double GetZ()
        {
            return z;
        }
    }
    class SignCurve
    {
        private const int NUM_RESAMPLE_POINTS = 20;
        private const int NUM_PARTS_R = 5;
        private const int NUM_PARTS_FI = 12;
        private const int NUM_PARTS_THETA = 6;
        private const double FI_0 = 2 * Math.PI / NUM_PARTS_FI;
        private const double THETA_0 = Math.PI / NUM_PARTS_THETA;
        private const int SC_FEATURE_DIM = NUM_PARTS_R * NUM_PARTS_FI * NUM_PARTS_THETA;
        private const double GAUSS_SIGMA = 1.0;
        private List<point> jointsCoor = new List<point>();
        private List<point> jointsCoor_resample = new List<point>();
        private List<point> jointsCoor_resample_norm = new List<point>();
        private List<double[]> scFeatures = new List<double[]>();
        private double curveLength;

        public SignCurve()
        {
            this.curveLength = 0;
            for (int i = 0; i < 10; i++)
            {
                Random ra = new Random();
                point p = new point(i, i, i);
                this.jointsCoor.Add(p);
            }
            CalculateLength();
        }

        public void WriteSCFeatures()
        {

        }

        public void GetShapeContext_3D()
        {
            Resample(NUM_RESAMPLE_POINTS);
            int numPoints = this.jointsCoor_resample.Count();
            double[] curve_x = new double[numPoints];
            double[] curve_y = new double[numPoints];
            double[] curve_z = new double[numPoints];
            for (int i = 0; i < numPoints; i++)
            {
                point p = this.jointsCoor_resample[i];
                curve_x[i] = p.GetX() / this.curveLength;
                curve_y[i] = p.GetY() / this.curveLength;
                curve_z[i] = p.GetZ() / this.curveLength;

                point p_norm = new point(curve_x[i], curve_y[i], curve_z[i]);
                this.jointsCoor_resample_norm.Add(p_norm);
            }

            double cMax_x = curve_x.Max();
            double cMax_y = curve_y.Max();
            double cMax_z = curve_z.Max();
            double cMin_x = curve_x.Min();
            double cMin_y = curve_y.Min();
            double cMin_z = curve_z.Min();
            double maxLen = Math.Sqrt(Math.Pow(cMax_x - cMin_x, 2) + Math.Pow(cMax_y - cMin_y, 2) + Math.Pow(cMax_z - cMin_z, 2));
            double R_normCoff = Math.Sqrt(Math.Pow(cMax_x - cMin_x, 2) + Math.Pow(cMax_y - cMin_y, 2) + Math.Pow(cMax_z - cMin_z, 2));
            double[] logLenR = new double[NUM_PARTS_R];

            for (int i = 0; i < NUM_PARTS_R; i++)
            {
                logLenR[i] = Math.Pow(0.5, NUM_PARTS_R - i - 1) * R_normCoff;
            }

            /*
            for (int i = 0; i < NUM_PARTS_R - 1; i++)
            {
                logLenR[i] = Math.Log((i + 1) * R_normCoff / NUM_PARTS_R);
            }
            logLenR[NUM_PARTS_R - 1] = double.PositiveInfinity;
            */
            for (int i = 0; i < numPoints; i++)
            {
                point p_ref = this.jointsCoor_resample_norm[i];
                List<point> tempCurve = new List<point>(this.jointsCoor_resample_norm.ToArray());
                tempCurve.RemoveAt(i);
                double[] feature_one_point = new double[SC_FEATURE_DIM];
                double normCoff = 0;
                for (int j = 0; j < tempCurve.Count(); j++)
                {
                    point p_dst = tempCurve[j];
                    int idx_bin = GetBinIdx(p_ref, p_dst, logLenR);
                    double dist = EuclideanDist(p_ref, p_dst);
                    double gauss_weight = GaussianFunc(GAUSS_SIGMA, 0, dist);
                    feature_one_point[idx_bin] += gauss_weight;
                    normCoff += gauss_weight;
                }
                for (int j = 0; j < SC_FEATURE_DIM; j++)
                {
                    feature_one_point[j] /= normCoff;
                }
                this.scFeatures.Add(feature_one_point);
            }
        }

        private double GaussianFunc(double sigma, double mu, double x)
        {
            double y = (1.0 / (2 * Math.PI * sigma)) * Math.Exp(-Math.Pow(x - mu, 2) / (2 * sigma * sigma));
            return y;
        }

        private int GetBinIdx(point p_ref, point p_dst, double[] logLenR)
        {
            double dx = p_dst.GetX() - p_ref.GetX();
            double dy = p_dst.GetY() - p_ref.GetY();
            double dz = p_dst.GetZ() - p_ref.GetZ();
            double dr = Math.Sqrt(dx * dx + dy * dy + dz * dz);

            double fi = (Math.Atan2(dy, dx) + 2 * Math.PI) % (2 * Math.PI);
            double theta = Math.Acos(dz / dr);

            int idx_r = 0;
            int idx_fi = (int)(fi / FI_0);
            int idx_theta = (int)(theta / THETA_0);

            for (int i = 0; i < NUM_PARTS_R; i++)
            {
                if (dr <= logLenR[i])
                {
                    idx_r = i;
                    break;
                }
            }

            int idx_bin = (idx_fi * NUM_PARTS_THETA + idx_theta) * NUM_PARTS_R + idx_r;
            return idx_bin;
        }

        public List<point> GetJointsCoor()
        {
            return this.jointsCoor;
        }

        public List<point> GetJointsCoor_resample()
        {
            return this.jointsCoor_resample;
        }

        public void PrintList(List<point> L)
        {
            foreach (point p in L)
            {
                Console.WriteLine("{0:F4} {1:F4} {2:F4}", p.GetX(), p.GetY(), p.GetZ());
            }
        }

        public List<point> Resample(int N)
        {
            List<point> rawCurve = new List<point>(this.jointsCoor.ToArray());
            int numPoints = rawCurve.Count();
            double I = this.curveLength / (N - 1);
            double D = 0;
            jointsCoor_resample.Add(rawCurve[0]);
            for (int i = 1; i < rawCurve.Count(); i++)
            {
                point p1 = rawCurve[i - 1];
                point p2 = rawCurve[i];
                double d = EuclideanDist(p1, p2);
                if ((D + d) >= I)
                {
                    double qx = p1.GetX() + ((I - D) / d) * (p2.GetX() - p1.GetX());
                    double qy = p1.GetY() + ((I - D) / d) * (p2.GetY() - p1.GetY());
                    double qz = p1.GetZ() + ((I - D) / d) * (p2.GetZ() - p1.GetZ());
                    point newPoint = new point(qx, qy, qz);
                    this.jointsCoor_resample.Add(newPoint);
                    rawCurve.Insert(i, newPoint);
                    D = 0;
                }
                else
                {
                    D = D + d;
                }
            }
            return this.jointsCoor_resample;
        }

        private double CalculateLength()
        {
            int numPoints = this.jointsCoor.Count();
            double length = 0;
            for (int i = 1; i < numPoints; i++)
            {
                point p1 = this.jointsCoor[i - 1];
                point p2 = this.jointsCoor[i];
                length += EuclideanDist(p1, p2);
            }
            curveLength = length;
            return length;
        }

        private double EuclideanDist(point p1, point p2)
        {
            double dx2 = Math.Pow(p2.GetX() - p1.GetX(), 2);
            double dy2 = Math.Pow(p2.GetY() - p1.GetY(), 2);
            double dz2 = Math.Pow(p2.GetZ() - p1.GetZ(), 2);
            double dist = Math.Sqrt(dx2 + dy2 + dz2);
            return dist;
        }
    }
}
