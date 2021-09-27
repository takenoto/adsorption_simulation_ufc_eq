library discretization.dart;

class UnidimensionalDiscretization{
  final List<double> values;
  UnidimensionalDiscretization(this.values);

   UnidimensionalDiscretization clone(){
     return UnidimensionalDiscretization(
       this.values.map((e) => e).toList()
     );
   }
}

class BidimensionalDiscretization{
  final List<List<double>> values;
  BidimensionalDiscretization(this.values);
  factory BidimensionalDiscretization.zeros(int d1, int d2){
    final v = List.generate(d1, (i) => List.generate(d2, (i) => 0.0));
    return BidimensionalDiscretization(v);
  }
}

class TridimensionalDiscretization{
  final List<List<List<double>>> values;
  TridimensionalDiscretization(this.values);
  factory TridimensionalDiscretization.zeros(int d1, int d2, int d3){
    final v = List.generate(d1, (i) => List.generate(d2, (i) => 
                                                    List.generate(d3, (i) => .0)));
    return TridimensionalDiscretization(v);
  }
}