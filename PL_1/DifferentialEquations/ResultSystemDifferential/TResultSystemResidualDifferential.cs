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
    public class TSystemResidualResultDifferential : TSystemResultDifferential
    {
        /// <summary>
        /// Точки решения системы дифф. уравнений
        /// </summary>
        public new List<TPointSystemResidualDifferential> SystemPoints { get; set; } = new List<TPointSystemResidualDifferential>();
        //-------------------------------------------------------
    }
}
