// Унифицированный результат решения системы дифференциальных уравнений с n переменными и n неизвестными
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
//*********************************************************
namespace StandartHelperLibrary.MathHelper
{
    /// <summary>
    /// Унифицированный результат решения системы дифференциальных уравнений с n переменными и n неизвестными
    /// </summary>
    public class TSystemResidualResultDifferential
    {
        /// <summary>
        /// Исходная система уравнений и настройки решателя
        /// </summary>
        public ISystemDifferentialEquation Equation { get; set; }
        /// <summary>
        /// Точки решения системы дифф. уравнений
        /// </summary>
        public List<TPointSystemResidualDifferential> SystemPoints { get; set; } = new List<TPointSystemResidualDifferential>();
        //-------------------------------------------------------
        /// <summary>
        /// В текстовое представление
        /// </summary>
        /// <returns>Текстовое представление</returns>
        public override string ToString()
        {
            string Out = "";
            foreach (var P in SystemPoints)
            {
                Out = Out + P.ToString() + " \r\n";
            }
            return Out;
        }
    }
}
