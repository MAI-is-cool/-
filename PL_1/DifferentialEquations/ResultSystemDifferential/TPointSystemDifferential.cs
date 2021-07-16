// Точки решения системы дифференциальных уравнений с n переменными и n неизвестными
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
//*********************************************************
namespace StandartHelperLibrary.MathHelper
{
    /// <summary>
    /// Точки решения системы дифференциальных уравнений с n переменными и n неизвестными
    /// </summary>
    /// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    /// ЗДЕСЬ НАДО ПРИДУМАТЬ ХОРОШИЙ ВЫВОД
    public class TPointSystemDifferential
    {
        /// <summary>
        /// Номер итерации
        /// </summary>
        public int IndexIteration { get; set; } = 0;
        /// <summary>
        /// Значения (Х) передаваемые в уравнение 
        /// </summary>
        public double X;
        /// <summary>
        /// Результат расчета в точке
        /// </summary>
        public double[] Result { get; set; }
        /// <summary>
        /// Коэффициенты решения
        /// </summary>
        public List<double[]> Coeffs { get; set; }
//-------------------------------------------------------
        /// <summary>
        /// В текстовое представление
        /// </summary>
        /// <returns>Текстовое представление</returns>
        public override string ToString()
        {
            string strResult = "Point №" + (IndexIteration+1).ToString();
            for (int i = 0; i < Coeffs.Count; i++)
            {
                strResult += "\nK" + (i + 1).ToString() + ":  ";
                for (int j = 0; j < Coeffs[i].Length; j++)
                    strResult += Coeffs[i][j].ToString() + "; ";
            }
            strResult += "\nX:   " + X.ToString() + "\n";
            for (int i = 0; i < Result.Length; i++)
                strResult += "Y" + (i+1).ToString() + ":  " + Result[i].ToString() + "\n";
            

            return strResult;
        }
        //-------------------------------------------------------
    }
}
