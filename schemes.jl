abstract type SplittingScheme end
  
struct MagnusScheme
    c::Array{Float64,2}
    a::Array{Float64,2}
    order::Integer
    name::AbstractString
end

struct EmbeddedScheme <: SplittingScheme
    scheme::Tuple
    scheme2::Tuple
    operator_sequence::String
    order::Integer
    name::AbstractString
end

struct MilneScheme <: SplittingScheme
    scheme::Tuple
    scheme2::Tuple
    operator_sequence::String
    order::Integer
    k::Float64
    name::AbstractString
end

struct PalindromicScheme <: SplittingScheme
    scheme::Tuple
    operator_sequence::String
    order::Integer
    name::AbstractString
end

struct DefectBasedScheme <: SplittingScheme
    scheme::Tuple
    operator_sequence::String
    order::Integer
    name::AbstractString
end

struct EmbeddedMagnusScheme
    scheme::MagnusScheme
    scheme2::MagnusScheme
    order::Integer
    name::AbstractString
end

struct ClassicMagnusScheme
    c::Array{Float64,1}
    order::Integer
    name::AbstractString
end

########################################################################
# _____________________EmbeddedScheme________________________________________________________________#
# Emb 3/2 AK S
const EMB32AKs = EmbeddedScheme(
          (0.268330095781759925,  0.919661523017399857,
           -0.187991618799159782, -0.187991618799159782,
            0.919661523017399857,  0.268330095781759925),

          (0.268330095781759925,  0.683368274569431595,
            0.731669904218240075,  0.316631725430568405),

            "AB", 2 ,"EMB 3/2 AK s")

# Emb 4/3 AK p
const EMB43AKp = EmbeddedScheme(
          (0.125962888700250514,  0.333588446797901933,
            0.751193431379145450, -0.338296598434303506,
            0.127551831557005609,  0.127551831557005609,
           -0.338296598434303506,  0.751193431379145450,
            0.333588446797901933,  0.125962888700250514),

          (0.125962888700250514,  0.333588446797901933,
            0.751193431379145450, -0.338296598434303506,
            0.0,  0.261153550449697153,
           -0.242703571757396124,  0.596114052266110425,
            0.365547251678000160,  0.147440548920593995),

            "AB", 3 ,"EMB 4/3 AK p")
            
# Emb 4/3 AK s
const EMB43AKs = EmbeddedScheme(
          (0.267171359000977615, -0.361837907604416033,
           -0.0338279096695056672,  0.861837907604416033,
            0.5333131013370561044,  0.861837907604416033,
           -0.0338279096695056672, -0.361837907604416033,
            0.267171359000977615,  0.0),

          (0.267171359000977615, -0.361837907604416033,
           -0.0338279096695056672,  0.861837907604416033,
            0.5333131013370561044,  0.395088376480991403,
            0.267171359000977615, -0.361837907604416033,
           -0.0338279096695056672,  0.466749531123424630),
           
           "AB", 3 ,"EMB 4/3 AK s")

# Emb 4/3 M/AK 
const EMB43MAK = EmbeddedScheme(
          (0.0935003487263305760,  0.439051727817158558,
           -0.0690943698810950380, -0.136536314071511211,
            0.475594021154764462,  0.394969172508705306,
            0.475594021154764462, -0.136536314071511211,
           -0.0690943698810950380,  0.439051727817158558,
            0.0935003487263305760,  0.0),

          (0.0935003487263305760,  0.439051727817158558,
           -0.0690943698810950380, -0.136536314071511211,
            0.475594021154764462,  0.394969172508705306,
            0.344129998965041583,  0.439051727817158558,
           -0.0690943698810950380, -0.136536314071511211,
            0.224964370916053455,  0.0),
           
            "AB", 3 ,"EMB 4/3 M/AK")

# Emb 4/3 BM PRK/A 
const EMB43BMPRKA = EmbeddedScheme(
          (0.0792036964311954608,  0.209515106613361891,
            0.353172906049773948, -0.143851773179818077,
           -0.0420650803577191948,  0.434336666566456186,
            0.219376955753499572,  0.434336666566456186,
           -0.0420650803577191948, -0.143851773179818077,
            0.353172906049773948,  0.209515106613361891,
            0.0792036964311954608,  0.0),

          (0.0792036964311954608,  0.209515106613361891,
            0.353172906049773948,  0.634279607840446390,
           -0.351497633364616918, -0.0576123055740887430,
            0.919121030883647509,  0.213817591120280462),
           
            "AB", 3 ,"EMB 4/3 BM PRK/A ")
# Emb 5/4 AK (ii) 
const EMB54AKii = EmbeddedScheme(
          (0.475018345144539497, -0.402020995028838599,
            0.021856594741098449,  0.345821780864741783,
           -0.334948298035883491,  0.400962967485371350,
            0.512638174652696736,  0.980926531879316517,
           -0.011978701020553904, -1.362064898669775624,
           -0.032120004263046859,  0.923805029000837468,
            0.369533888781149572,  0.112569584468347105,),

          (0.475018345144539497, -0.402020995028838599,
            0.021856594741098449,  0.345821780864741783,
           -0.334948298035883491,  0.400962967485371350,
            0.512638174652696736,  0.980926531879316517,
           -0.00596874002994121298, -1.280062271930040874,
           -0.0418443254169191836,  0.709119365029440952,
            0.0616955797702736790,  0.132683037231661766,
            0.31155266917413552658,  0.112569584468347105),
           
            "AB", 4 ,"EMB 5/4 AK (ii)")
            
# PP 3/4 A
const PP34A = PalindromicScheme( 
          (0.268330095781759925,  0.919661523017399857, 
           -0.187991618799159782, -0.187991618799159782, 
            0.919661523017399857,  0.268330095781759925),

            "AB", 3 ,"PP 3/4 A")
            
# PP 5/6 A
const PP56A = PalindromicScheme(
          (0.201651044312324230,  0.578800656272664932, 
            0.562615975356569200,  0.273128836056524479, 
            0.253874038247554845, -0.102733803148432142, 
           -0.835351693190370636,  0.068014946093165092, 
            0.068014946093165092, -0.835351693190370636,
           -0.102733803148432142,  0.253874038247554845, 
            0.273128836056524479,  0.562615975356569200, 
            0.578800656272664932,  0.201651044312324230),
            
            "AB", 5 ,"PP 5/6 A")
            
# Symm-Milne-32
const SymmMilne32 = MilneScheme(
          (0.190983005625052576,  0.500000000000000000,
            0.618033988749894848,  0.500000000000000000,
            0.190983005625052576,  0.0),

          (0.292893218813452475,  0.707106781186547525,
            0.707106781186547525,  0.292893218813452475),
            
            "AB", 2 ,1.5081045295901757 ,"Symm-Milne-32")

# _____________________defectbased________________________________________________________________#
const Strang = DefectBasedScheme(
            (0.5, 1.0,
             0.5, 0.0),"AB",  2 ,"Strang")
const A106 = DefectBasedScheme(
            (0.0951762545417740527, 0.666296893997707801,
            -0.127950285523686779,  0.0246189009521050871,
             0.105972953453251131, -0.410725533617951132,
             0.4482222766008274842, 0.657729262050913178,
            -0.0214211990721658889, -0.8758390467655498682,
            -0.0214211990721658889, 0.657729262050913178,
             0.4482222766008274842, -0.410725533617951132,
             0.105972953453251131,  0.0246189009521050871,
            -0.127950285523686779,  0.666296893997707801,
             0.0951762545417740527, 0.000000000000000000),
             "AB", 6 ,"A 10-6")
const BM116PRK = DefectBasedScheme(
             (0.0502627644003922, 0.148816447901042,
              0.413514300428344, -0.132385865767784,
              0.0450798897943977, 0.067307604692185,
             -0.188054853819569,  0.432666402578175,
              0.541960678450780, -0.016404589403618,
             -0.7255255585086898, -0.016404589403618,
              0.541960678450780,  0.432666402578175,
             -0.188054853819569,  0.067307604692185,
              0.0450798897943977, -0.132385865767784,
              0.413514300428344,  0.148816447901042,
              0.0502627644003922, 0.000000000000000),
              "AB", 6 ,"BM 11-6 PRK")
             
const EMB43AKp_D4 = DefectBasedScheme(EMB43AKp.scheme, "AB", EMB43AKp.order + 1, EMB43AKp.name * " (defect 4)")
const EMB43AKp_D3 = DefectBasedScheme(EMB43AKp.scheme2, "AB", EMB43AKp.order,  EMB43AKp.name * " (defect 3)")
const EMB43AKs_D4 = DefectBasedScheme(EMB43AKs.scheme, "AB", EMB43AKs.order + 1, EMB43AKs.name * " (defect 4)")
const EMB43AKs_D3 = DefectBasedScheme(EMB43AKs.scheme2, "AB", EMB43AKs.order,  EMB43AKs.name * " (defect 3)")
const EMB43MAK_D4 = DefectBasedScheme(EMB43MAK.scheme, "AB", EMB43MAK.order + 1, EMB43MAK.name * " (defect 4)")
const EMB43MAK_D3 = DefectBasedScheme(EMB43MAK.scheme2, "AB", EMB43MAK.order,  EMB43MAK.name * " (defect 3)")
const EMB43BMPRKA_D4 = DefectBasedScheme(EMB43BMPRKA.scheme, "AB", EMB43BMPRKA.order + 1, EMB43BMPRKA.name * " (defect 4)")
const EMB43BMPRKA_D3 = DefectBasedScheme(EMB43BMPRKA.scheme2, "AB", EMB43BMPRKA.order, EMB43BMPRKA.name * " (defect 3)")
const EMB54AKii_D5 = DefectBasedScheme(EMB54AKii.scheme, "AB", EMB54AKii.order + 1, EMB54AKii.name * " (defect 5)")
const EMB54AKii_D4 = DefectBasedScheme(EMB54AKii.scheme2, "AB", EMB54AKii.order, EMB54AKii.name * " (defect 4)")

const PP34A_D = DefectBasedScheme(PP34A.scheme, "AB", PP34A.order, PP34A.name * " (defect)")
const PP56A_D = DefectBasedScheme(PP56A.scheme, "AB", PP56A.order, PP56A.name * " (defect)")

const SymmMilne32_D = DefectBasedScheme(SymmMilne32.scheme, "AB", SymmMilne32.order, SymmMilne32.name * " (defect)")
# _____________________Commutator Free MagnusSchemes________________________________________________________________#
const Magnus2 = MagnusScheme(
            reshape([0.5], 1, 1),
            reshape([1.0], 1, 1),
            2 , "Magnus2")
const Magnus4 = MagnusScheme(
            repeat(transpose([0.5 - sqrt(3.0) / 6, 0.5 + sqrt(3.0) / 6]), 2),
            [0.25 + sqrt(3.0) / 6  0.25 - sqrt(3.0) / 6;
             0.25 - sqrt(3.0) / 6  0.25 + sqrt(3.0) / 6;],
            4 , "Magnus4")
const Magnus4O = MagnusScheme(
            repeat(transpose([1 / 2 - sqrt(3 / 20), 1 / 2, 1 / 2 + sqrt(3 / 20)]), 3),
            [ 37 / 240 + 10 / 87 * sqrt(5 / 3) -1 / 30  37 / 240 - 10 / 87 * sqrt(5 / 3);
             -11 / 360                 23 / 45 -11 / 360;
              37 / 240 - 10 / 87 * sqrt(5 / 3) -1 / 30  37 / 240 + 10 / 87 * sqrt(5 / 3);],
            4 , "Magnus4Optimized")
const Magnus4Brunner = MagnusScheme(
            [ 0.13471638221006548022e1   0.19019724244985505314e0;
              0.80980275755014494686e0  -0.34716382210065480219e0;],
            [-0.10169081889249764835e-1  0.51016908188924976483e0;
              0.51016908188924976483e0  -0.10169081889249764835e-1;],
            4, "Magnus4Brunner")
const Magnus6_4 = MagnusScheme(
            repeat(transpose([1 / 2 - sqrt(15) / 10, 1 / 2, 1 / 2 + sqrt(15) / 10]), 4),
            [-0.20052856057448226890e0  0.18713900774428756529e1 -0.59100909048596250154e0; 
              0.32223594293373734800e0 -0.16491678552206534307e1  0.74707948590448520023e0;
              0.74707948590448520023e0 -0.16491678552206534307e1  0.32223594293373734800e0; 
             -0.59100909048596250154e0  0.18713900774428756529e1 -0.20052856057448226890e0;],
            6 , "Magnus6_4")
const Magnus6_6 = MagnusScheme(
            repeat(transpose([1 / 2 - sqrt(15) / 10, 1 / 2, 1 / 2 + sqrt(15) / 10]), 6),
            [0.21583899697576773803e0  -0.76717964591551449767e-1  0.20878967615783711741e-1;
            -0.80897796320852999368e-1 -0.17874721753715766251e0   0.32263366431047360188e-1;
             0.18062846005583006946e0   0.47768740435093133450e0  -0.90934216979798102268e-1;
            -0.90934216979798102268e-1  0.47768740435093133450e0   0.18062846005583006946e0;
             0.32263366431047360188e-1 -0.17874721753715766251e0  -0.80897796320852999368e-1;
             0.20878967615783711741e-1 -0.76717964591551449767e-1  0.21583899697576773803e0;],
            6 , "Magnus6_6")
const MagnusB3 = MagnusScheme(
            [0.0 0.0; 0.710937500000                   -0.85],
            [1 / 3 0.0; 0.683350016683350016683350016683 -0.0166833500166833500166833500167],
            3 , "MagnusB3")
const MagnusB2 = MagnusScheme(
            [0.0 0.0; 3 / 4 -0.85],
            [1 / 3 0.0; 2 / 3  0.0],
            2 , "MagnusB2")
const MagnusB3_new = MagnusScheme(
            [0.19945632824746859504924911091   0.19945632824746859504924911091;
             0.787136184478554178003475703186  0.122999272714433843887903890536],
            [0.554550577274841171530710728519  0.0;
             0.503813923117446911959772513908 -0.0583645003922880834904832424241],
            3 , "MagnusB3_new")
const MagnusB2_new = MagnusScheme(
            [0.19945632824746859504924911091  0.19945632824746859504924911091;
             0.0                              0.799380164239681074819809203679],
            [0.554550577274841171530710728519 0.0;
            -0.041667215200939925165248146774 0.487116637926098753634537418255],
            2 , "MagnusB2_new")
const MagnusB4 = MagnusScheme(
            reshape([ 0.675603595979828817023843904485;
              0.500000000000000000000000000000;
              0.324396404020171182976156095515 ],3,1),
            reshape([ 1.35120719195965763404768780897;
             -1.70241438391931526809537561794;
              1.35120719195965763404768780897],3,1),
            4 , "MagnusB4")
const MagnusEMB32 = EmbeddedMagnusScheme(MagnusB3, MagnusB2, 3, "Magnus EMB 3/2")
const MagnusEMB32_new = EmbeddedMagnusScheme(MagnusB3_new, MagnusB2_new, 3, "Magnus EMB 3/2_new")

const CF4oH = MagnusScheme(
# optimized with respect to LEM with equal weight:
# [ 0.311314623379755386999882845054  -0.027985859584027834823234100810  0.007787203484903714984134658507
# -0.041324049086881324206239725785   0.500416163612500114090912646066 -0.041324049086881324206239725788
#  0.0077872034849037149841346585064 -0.027985859584027834823234100810  0.311314623379755386999882845055],    
#
# optimized with respect to scaled LEM with weights 
# [A1,[A1,[A1,A2]]]: 1, [A2,[A1,A2]]: 0.031, [A1,[A1,A3]]: 0.077, [A2,A3]: 0.00013
   repeat(transpose([1 / 2 - sqrt(15) / 10, 1 / 2, 1 / 2 + sqrt(15) / 10]), 3),
   [  0.302146842308616954258187683416  -0.030742768872036394116279742324  0.004851603407498684079562131338;
  -0.029220667938337860559972036973   0.505929982188517232677003929089 -0.029220667938337860559972036973;
   0.004851603407498684079562131337  -0.030742768872036394116279742324  0.302146842308616954258187683417],
     4,"CF4oH")

const CF6n = MagnusScheme(
      repeat(transpose([1 / 2 - sqrt(15) / 10, 1 / 2, 1 / 2 + sqrt(15) / 10]), 4),
  [ 7.9124225942889763e-01 -8.0400755305553218e-02  1.2765293626634554e-02;
   -4.8931475164583259e-01  5.4170980027798808e-02 -1.2069823881924156e-02;
   -2.9025638294289255e-02  5.0138457552775674e-01 -2.5145341733509552e-02;
    4.8759082890019896e-03 -3.0710355805557892e-02  3.0222764976657693e-01],
    
  6,"CF6n")
# _____________________ClassicMagnusSchemes_______________________________________________#
ClassicMagnus4 = ClassicMagnusScheme(
            [0.5 - sqrt(3) / 6; 0.5 + sqrt(3) / 6],
            4 , "ClassicMagnus4")
# _____________________ABC________________________________________________________________#
const LieTrotterABC = DefectBasedScheme(
            (1.0, 1.0, 1.0), "ABC", 2 ,"LieTrotter ABC")
            
const StrangABC = DefectBasedScheme(
            (0.0, 0.0, 0.5,
             0.0, 0.5, 0.0,
             1.0, 0.5, 0.5), "ABC", 2 ,"Strang ABC")
             
w1 = 1 / (2 - 2^(1 / 3))
w0 = 1 - 2 * w1
const Y74ABC = DefectBasedScheme(
            (0.0, 0.0, w1 / 2,
             0.0, w1 / 2, 0.0,
             w1, w1 / 2, (w1 + w0) / 2,
             0.0, w0 / 2, 0.0,
             w0, w0 / 2, (w1 + w0) / 2,
             0.0, w1 / 2, 0.0,
             w1, w1 / 2, w1 / 2), 
             "ABC", 4 , "Y 7-4")
w1 = -1.177679984178871007 
w2 = 0.235573213359358134 
w3 = 0.784513610477557264 
w0 = 1 - 2 * (w1 + w2 + w3)
const AY156ABC = DefectBasedScheme(
            (0.0, 0.0, w3 / 2,
             0.0, w3 / 2, 0.0,
             w3, w3 / 2, (w2 + w3) / 2,
             0.0, w2 / 2, 0.0,
             w2, w2 / 2, (w1 + w2) / 2,
             0.0, w1 / 2, 0.0,
             w1, w1 / 2, (w0 + w1) / 2,
             0.0, w0 / 2, 0.0,
             w0, w0 / 2, (w0 + w1) / 2,
             0.0, w1 / 2, 0.0,
             w1, w1 / 2, (w1 + w2) / 2,
             0.0, w2 / 2, 0.0,
             w2, w2 / 2, (w2 + w3) / 2,
             0.0, w3 / 2, 0.0,
             w3, w3 / 2, w3 / 2),
             "ABC", 6 , "AY 15-6")
