#include "csvIO.hpp"
#include "file_explorer.hpp"
#include "linalg.hpp"
#include "logical_operators.hpp"
#include "matrix.hpp"
#include "mmath.hpp"
#include "optimization.hpp"
#include "polynomial.hpp"
#include "regression.hpp"
#include "SDL_interface.hpp"
#include "statistics.hpp"
#include "test_functions.hpp"

#include "MachineLearning/nnmf.hpp"
#include "MachineLearning/pca.hpp"
#include "MachineLearning/scaler.hpp"
#include "MachineLearning/Clustering/affinitypropagation.hpp"
#include "MachineLearning/Clustering/gmm.hpp"
#include "MachineLearning/Clustering/hierarchicalclustering.hpp"
#include "MachineLearning/Clustering/kmeans.hpp"
#include "MachineLearning/Clustering/kmedoids.hpp"
#include "MachineLearning/Clustering/meanshift.hpp"
#include "MachineLearning/SupervisedLearning/bernoulli_naive_bayes.hpp"
#include "MachineLearning/SupervisedLearning/dtree.hpp"
#include "MachineLearning/SupervisedLearning/gaussian_naive_bayes.hpp"
#include "MachineLearning/SupervisedLearning/knn.hpp"
#include "MachineLearning/SupervisedLearning/lda.hpp"
#include "MachineLearning/SupervisedLearning/mlp.hpp"
#include "MachineLearning/SupervisedLearning/perceptron.hpp"
#include "MachineLearning/SupervisedLearning/qda.hpp"

#include "ImageProcessing/color.hpp"
#include "ImageProcessing/fourier.hpp"
#include "ImageProcessing/hough.hpp"
#include "ImageProcessing/imgconverter.hpp"
#include "ImageProcessing/linear_filtering.hpp"
#include "ImageProcessing/morphology.hpp"
#include "ImageProcessing/nonlinear_filtering.hpp"
#include "ImageProcessing/segmentation.hpp"
#include "ImageProcessing/shape_features.hpp"
#include "ImageProcessing/snake.hpp"

#include "ImageProcessing/ImageIO/imageIO.hpp"
#include "ImageProcessing/ImageIO/bmpIO.hpp"
#include "ImageProcessing/ImageIO/gifIO.hpp"
#include "ImageProcessing/ImageIO/icoIO.hpp"
#include "ImageProcessing/ImageIO/jpegIO.hpp"
#include "ImageProcessing/ImageIO/pngIO.hpp"
#include "ImageProcessing/ImageIO/pnmIO.hpp"
#include "ImageProcessing/ImageIO/tgaIO.hpp"
#include "ImageProcessing/ImageIO/tiffIO.hpp"

#include <complex>
void test_clustering()
{
     ///////////test clustering
    //create data
    std::size_t n_samples = 240;
    std::size_t n_features = 2;
    std::size_t n_classes = 4;
    Matrix<float> Data(n_samples, n_features);
    Matrix<std::size_t> Labels(n_samples, 1);
    Matrix<float> means = {{10.0f, 20.0f, 30.0f, 40.0f},
                            {15.0f, 40.0f, 5.0f, 30.0f}};
    Matrix<float> stddevs = ones<float>(n_features, n_classes);
    Matrix<float> results;

    for(std::size_t i=0;i<n_classes;++i)
    {
        for(std::size_t j=0;j<n_features;++j)
            Data.setSubmat(i*n_samples/n_classes, j, randn(n_samples/n_classes, 1, means(j, i), stddevs(j, i)));
        Labels.setSubmat(i*n_samples/n_classes, 0, full<std::size_t>(n_samples/n_classes, 1, i));
    }

    // test Affinity Propagation
    std::cout<<"Affinity Propagation"<<std::endl;
    AffinityPropagation clf_affpro;
    results = Matrix<float>(clf_affpro.fit_predict(Data));
    results.reshape(n_classes, n_samples/n_classes);
    std::cout<<(histogram(Matrix<std::size_t>(results))).getCols(0, n_classes)<<std::endl;
    std::cout<<axismean(results, 1)<<std::endl;
    std::cout<<axisstdev(results, 1)<<std::endl;
    std::cout<<clf_affpro.getCenters()<<std::endl;

    // test Kmeans
    std::cout<<"Kmeans"<<std::endl;
    Kmeans clf(n_classes);
    results = Matrix<float>(clf.fit_predict(Data));
    results.reshape(n_classes, n_samples/n_classes);
    std::cout<<(histogram(Matrix<std::size_t>(results))).getCols(0, n_classes)<<std::endl;
    std::cout<<axismean(results, 1)<<std::endl;
    std::cout<<axisstdev(results, 1)<<std::endl;
    std::cout<<clf.getCenters()<<std::endl;

    // test, Kmedoids
    std::cout<<"Kmedoids"<<std::endl;
    Kmedoids kmed_clf(n_classes);
    results = Matrix<float>(kmed_clf.fit_predict(Data));
    results.reshape(n_classes, n_samples/n_classes);
    std::cout<<(histogram(Matrix<std::size_t>(results))).getCols(0, n_classes)<<std::endl;
    std::cout<<axismean(results, 1)<<std::endl;
    std::cout<<axisstdev(results, 1)<<std::endl;
    std::cout<<kmed_clf.getCenters()<<std::endl;

    // test hierarchical clustering
    std::cout<<"hierarchical clustering"<<std::endl;
    HierarchicalClustering hier_clf(n_classes);
    results = Matrix<float>(hier_clf.fit_predict(Data));
    results.reshape(n_classes, n_samples/n_classes);
    std::cout<<(histogram(Matrix<std::size_t>(results))).getCols(0, n_classes)<<std::endl;
    std::cout<<axismean(results, 1)<<std::endl;
    std::cout<<axisstdev(results, 1)<<std::endl;


    //test meanshift
    std::cout<<"Mean Shift"<<std::endl;
    MeanShift clf_meanshift = MeanShift();
    results = Matrix<float>(clf_meanshift.fit_predict(Data));
    results.reshape(n_classes, n_samples/n_classes);
    std::cout<<(histogram(Matrix<std::size_t>(results))).getCols(0, n_classes)<<std::endl;
    std::cout<<axismean(results, 1)<<std::endl;
    std::cout<<axisstdev(results, 1)<<std::endl;

    //test gmm
    std::cout<<"diagonal GMM"<<std::endl;
    stddevs = {{1.0f, 3.0f, 1.0f, 1.0f},
                {2.0f, 1.5f, 2.5f, 2.0f}};
    for(std::size_t i=0;i<n_classes;++i)
    {
        for(std::size_t j=0;j<n_features;++j)
            Data.setSubmat(i*n_samples/n_classes, j, randn(n_samples/n_classes, 1, means(j, i), stddevs(j, i)));
        Labels.setSubmat(i*n_samples/n_classes, 0, full<std::size_t>(n_samples/n_classes, 1, i));
    }
    GMM clf_gmm(n_classes, "diag");
    clf_gmm.fit(Data);
    results = Matrix<float>(clf_gmm.fit_predict(Data));
    results.reshape(n_classes, n_samples/n_classes);
    std::cout<<(histogram(Matrix<std::size_t>(results))).getCols(0, n_classes)<<std::endl;
    std::cout<<axismean(results, 1)<<std::endl;
    std::cout<<axisstdev(results, 1)<<std::endl;
    std::cout<<clf_gmm.getCenters()<<std::endl;
    for(std::size_t i=0;i<n_classes;++i)
        std::cout<<clf_gmm.getCovariances()[i]<<std::endl;
    std::cout<<clf_gmm.getWeights()<<std::endl;

}

void test_supervisedlearning()
{
    //create data
    std::size_t n_samples = 2000;
    std::size_t n_features = 2;
    std::size_t n_classes = 4;
    Matrix<float> Data(n_samples, n_features);
    Matrix<std::size_t> Labels(n_samples, 1);
    Matrix<float> means = {{10.0f, 20.0f, 30.0f, 40.0f},
                            {15.0f, 40.0f, 5.0f, 30.0f}};

    Matrix<float> stddevs = ones<float>(n_features, n_classes);

    Matrix<float> results;

    for(std::size_t i=0;i<n_classes;++i)
    {
        for(std::size_t j=0;j<n_features;++j)
            Data.setSubmat(i*n_samples/n_classes, j, randn(n_samples/n_classes, 1, means(j, i), stddevs(j, i)));
        Labels.setSubmat(i*n_samples/n_classes, 0, full<std::size_t>(n_samples/n_classes, 1, i));
    }

    ///// test Dtree
    std::cout<<"Dtree"<<std::endl;
    Dtree tree;
    results = Matrix<float>(tree.fit_predict(Data, Labels));
    std::cout<<count_nonzero(Matrix<std::size_t>(results) == Labels)<<std::endl;
    tree.export_graphviz("tree.dot");

    ///// test k-NN
    std::cout<<"KNN"<<std::endl;
    KNN clf_knn(5);
    clf_knn.fit(Data, Labels);
    results = Matrix<float>(clf_knn.fit_predict(Data, Labels));
    std::cout<<count_nonzero(Matrix<std::size_t>(results) == Labels)<<std::endl;

    ///// test MLP
    std::cout<<"MLP"<<std::endl;
    StandardScaler ssc(1,1, 2);
    Data = ssc.fit_transform(Data);
    std::vector<std::size_t> param = {};
    MLP clf_mlp(param);
    results = Matrix<float>(clf_mlp.fit_predict(Data, Labels));
    clf_mlp.export_graphviz("test.dot");
    std::cout<<count_nonzero(Matrix<std::size_t>(results) == Labels)<<std::endl;

    ///// test gaussian_naives_bayes
    std::cout<<"Gaussian naives bayes"<<std::endl;
    GaussianNaiveBayes clf_gnb;
    clf_gnb.fit(Data, Labels);
    results = Matrix<float>(clf_gnb.fit_predict(Data, Labels));
    std::cout<<count_nonzero(Matrix<std::size_t>(results) == Labels)<<std::endl;
    std::cout<<clf_gnb.getCenters()<<std::endl;
    std::cout<<clf_gnb.getStddevs()<<std::endl;
    std::cout<<clf_gnb.getPriors()<<std::endl;

    ///// test LDA
    std::cout<<"LDA"<<std::endl;
    LDA clf_lda;
    clf_lda.fit(Data, Labels);
    results = Matrix<float>(clf_lda.fit_predict(Data, Labels));
    std::cout<<count_nonzero(Matrix<std::size_t>(results) == Labels)<<std::endl;
    std::cout<<clf_lda.getCenters()<<std::endl;
    std::cout<<clf_lda.getCovariances()<<std::endl;
    std::cout<<clf_lda.getPriors()<<std::endl;

    ///// test QDA
    std::cout<<"QDA"<<std::endl;
    QDA clf_qda;
    clf_qda.fit(Data, Labels);
    results = Matrix<float>(clf_qda.fit_predict(Data, Labels));
    std::cout<<count_nonzero(Matrix<std::size_t>(results) == Labels)<<std::endl;
    std::cout<<clf_qda.getCenters()<<std::endl;
    for(std::size_t i=0;i<n_classes;++i)
        std::cout<<clf_qda.getCovariances()[i]<<std::endl;
    std::cout<<clf_qda.getPriors()<<std::endl;
}

void test_morpho_binaire()
{
    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window* pWindow = SDL_CreateWindow("Image Processing Library",
                                            SDL_WINDOWPOS_UNDEFINED,
                                            SDL_WINDOWPOS_UNDEFINED,
                                            640,
                                            640,
                                            SDL_WINDOW_SHOWN);

    SDL_Renderer* renderer = SDL_CreateRenderer(pWindow, -1, SDL_RENDERER_ACCELERATED);

    for(int i=1;i<5;++i)
        show_matrix(renderer, circle(i));
    for(int i=1;i<5;++i)
        show_matrix(renderer, cross(i));
    for(int i=1;i<5;++i)
    {
        for(int j=1;j<5;++j)
            show_matrix(renderer, ellipse(i,j));
    }
    for(int i=1;i<5;++i)
        show_matrix(renderer, diamond(i));
    for(int i=1;i<5;++i)
        show_matrix(renderer, octagon(i));
    for(int i=1;i<5;++i)
    {
        for(int j=1;j<5;++j)
            show_matrix(renderer, rect(i,j));
    }
    for(int i=1;i<5;++i)
        show_matrix(renderer, square(i));

    /*Matrix<bool> conservative_smoothing(const Matrix<bool>&);
    Matrix<bool> contrast_enhancement(const Matrix<bool>&);
    Matrix<bool> HitOrMiss(const Matrix<bool>&, const Mask&, const Mask&);
    Matrix<bool> inner_gradient(const Matrix<bool>&, const Mask& = square(3));
    Matrix<bool> reconstruct(const Matrix<bool>&, const Matrix<bool>& , const Mask& = square(3));
    Matrix<bool> outer_gradient(const Matrix<bool>&, const Mask& = square(3));

    Matrix<std::size_t> label(const Matrix<bool>&, const Mask& = square(3));*/


    Matrix<bool> gray_img = {{0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {0,0,0,0,0,1,0,0,0,0,0,0,0,0},
                {0,0,0,0,0,1,1,1,1,0,0,0,0,0},
                {0,0,0,0,0,1,1,1,1,0,1,0,0,0},
                {0,0,0,0,0,0,1,1,1,1,0,0,0,0},
                {0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                {0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    };

    SDL_SetWindowTitle(pWindow, "original");
    show_matrix(renderer, gray_img);
    SDL_SetWindowTitle(pWindow, "erode");
    show_matrix(renderer, erode(gray_img));
    SDL_SetWindowTitle(pWindow, "dilate");
    show_matrix(renderer, dilate(gray_img));
    SDL_SetWindowTitle(pWindow, "open");
    show_matrix(renderer, open(gray_img));
    SDL_SetWindowTitle(pWindow, "close");
    show_matrix(renderer, close(gray_img));
    SDL_SetWindowTitle(pWindow, "BTH");
    show_matrix(renderer, BTH(gray_img));
    SDL_SetWindowTitle(pWindow, "WTH");
    show_matrix(renderer, WTH(gray_img));
    SDL_SetWindowTitle(pWindow, "gradient");
    show_matrix(renderer, gradient(gray_img));
    SDL_SetWindowTitle(pWindow, "inner gradient");
    show_matrix(renderer, inner_gradient(gray_img));
    SDL_SetWindowTitle(pWindow, "outer gradient");
    show_matrix(renderer, outer_gradient(gray_img));
    SDL_SetWindowTitle(pWindow, "convex hull");
    show_matrix(renderer, convex_hull(gray_img));
    SDL_SetWindowTitle(pWindow, "skeleton");
    show_matrix(renderer, skeleton(gray_img));
    SDL_SetWindowTitle(pWindow, "conservative smoothing");
    show_matrix(renderer, conservative_smoothing(gray_img));
    SDL_SetWindowTitle(pWindow, "contrast enhancement");
    show_matrix(renderer, contrast_enhancement(gray_img));
    SDL_SetWindowTitle(pWindow, "opening by reconstruction");
    gray_img = opening_by_reconstruction(gray_img);
    show_matrix(renderer, gray_img);
    gray_img = closing_by_reconstruction(gray_img);
    SDL_SetWindowTitle(pWindow, "closing by reconstruction");
    show_matrix(renderer, gray_img);

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(pWindow);
    SDL_Quit();
}

void test_morpho_gray()
{
    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window* pWindow = SDL_CreateWindow("Image Processing Library",
                                            SDL_WINDOWPOS_UNDEFINED,
                                            SDL_WINDOWPOS_UNDEFINED,
                                            640,
                                            480,
                                            SDL_WINDOW_SHOWN);

    SDL_Renderer* renderer = SDL_CreateRenderer(pWindow, -1, SDL_RENDERER_ACCELERATED);

    // test filtres non linéaires
    Matrix<Color> img = read_img("Images/Pixmap/clown_P5.pgm");
    Matrix<unsigned char> gray_img = color2grayimage(img);

    SDL_SetWindowTitle(pWindow, "original");
    show_matrix(renderer, gray_img);
    SDL_SetWindowTitle(pWindow, "erode");
    show_matrix(renderer, erode(gray_img));
    SDL_SetWindowTitle(pWindow, "dilate");
    show_matrix(renderer, dilate(gray_img));
    SDL_SetWindowTitle(pWindow, "open");
    show_matrix(renderer, open(gray_img));
    SDL_SetWindowTitle(pWindow, "close");
    show_matrix(renderer, close(gray_img));
    SDL_SetWindowTitle(pWindow, "BTH");
    show_matrix(renderer, BTH(gray_img));
    SDL_SetWindowTitle(pWindow, "WTH");
    show_matrix(renderer, WTH(gray_img));
    SDL_SetWindowTitle(pWindow, "hmin");
    show_matrix(renderer, hmin(gray_img));
    SDL_SetWindowTitle(pWindow, "hmax");
    show_matrix(renderer, hmax(gray_img));
    SDL_SetWindowTitle(pWindow, "median");
    show_matrix(renderer, median_filter(gray_img));

    /*SDL_SetWindowTitle(pWindow, "gradient");
    show_matrix(renderer, gradient(gray_img));*/
    SDL_SetWindowTitle(pWindow, "inner gradient");
    show_matrix(renderer, inner_gradient(gray_img));
    SDL_SetWindowTitle(pWindow, "outer gradient");
    show_matrix(renderer, outer_gradient(gray_img));
    /*SDL_SetWindowTitle(pWindow, "convex hull");
    show_matrix(renderer, convex_hull(gray_img));*/
    SDL_SetWindowTitle(pWindow, "skeleton");
    show_matrix(renderer, skeleton(gray_img));
    SDL_SetWindowTitle(pWindow, "conservative smoothing");
    show_matrix(renderer, conservative_smoothing(gray_img));
    SDL_SetWindowTitle(pWindow, "contrast enhancement");
    show_matrix(renderer, contrast_enhancement(gray_img));
    SDL_SetWindowTitle(pWindow, "opening by reconstruction");
    gray_img = opening_by_reconstruction(gray_img);
    show_matrix(renderer, gray_img);
    gray_img = closing_by_reconstruction(gray_img);
    SDL_SetWindowTitle(pWindow, "closing by reconstruction");
    show_matrix(renderer, gray_img);

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(pWindow);
    SDL_Quit();
}

void test_lecture_img()
{
    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window* pWindow = SDL_CreateWindow("Image Processing Library",
                                            SDL_WINDOWPOS_UNDEFINED,
                                            SDL_WINDOWPOS_UNDEFINED,
                                            640,
                                            480,
                                            SDL_WINDOW_SHOWN);

    SDL_Renderer* renderer = SDL_CreateRenderer(pWindow, -1, SDL_RENDERER_ACCELERATED);

    std::string extensions[] = { "bmp", "tga", "ico", "pbm", "pcx", "pgm", "ppm", "jpg", "gif", "png", "tiff"};

    std::for_each(
                  extensions,
                  extensions+9,
                  [renderer,pWindow](std::string s){auto cur_list = get_files_recursively(".", s);
                                                std::for_each(cur_list.begin(),
                                                  cur_list.end(),
                                                  [renderer,pWindow](std::string ss){std::cout<<ss<<std::endl;
                                                                            Matrix<Color> img = read_img(ss);
                                                                            SDL_SetWindowTitle(pWindow, ss.c_str());
                                                                            show_matrix(renderer, img);}
                                                  );
                                            }
                  );

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(pWindow);
    SDL_Quit();
}

void test_ecriture_img()
{
    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window* pWindow = SDL_CreateWindow("Image Processing Library",
                                            SDL_WINDOWPOS_UNDEFINED,
                                            SDL_WINDOWPOS_UNDEFINED,
                                            640,
                                            480,
                                            SDL_WINDOW_SHOWN);

    SDL_Renderer* renderer = SDL_CreateRenderer(pWindow, -1, SDL_RENDERER_ACCELERATED);



    //Test writer
    //Matrix<Color> img = read_img("Images/Detection_de_forme/bild4.png");
    Matrix<Color> img = read_img("Images/marinbas.ico");
    //img = read_img("Images/test_icone.ico");
    show_matrix(renderer, img);
    std::vector<std::string> extensions = { "bmp", "tga", "ico", "png", "pbm","pgm", "ppm", "jpg", "gif", "pcx", "tiff"};
    std::for_each(extensions.cbegin(), extensions.cend(), [img](const std::string& s){save_img("testwriter."+s, img);});
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(pWindow);
    SDL_Quit();
}

void test_fourier()
{
    // test FFT
    Matrix<float> data = sin(2.0f*float(PI)*0.2f*arange<float>(0.0f,  25.0f));
    std::cout<<data<<std::endl;
    std::cout<<FFT(data)<<std::endl;
    std::cout<<FFTshift(FFT(data))<<std::endl;
    std::cout<<iFFT(FFT(data))<<std::endl;
    std::cout<<FFTfreq(data.size())<<std::endl;

    Matrix<unsigned char> img = read_pgm("Images/Pixmap/antifeep_P2.pgm");

    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window* pWindow = SDL_CreateWindow("Image Processing Library",
                                            SDL_WINDOWPOS_UNDEFINED,
                                            SDL_WINDOWPOS_UNDEFINED,
                                            640,
                                            480,
                                            SDL_WINDOW_SHOWN);

    SDL_Renderer* renderer = SDL_CreateRenderer(pWindow, -1, SDL_RENDERER_ACCELERATED);
    show_matrix(renderer, img);
    Matrix<float> img_f(img);
    Matrix<std::complex<float> > fft_img = FFT(img_f);
    Matrix<std::complex<float> > ifft_img = iFFT(fft_img);

    std::cout<<FFTshift(fft_img)<<std::endl;
    Matrix<float> resultfft(apply<std::complex<float>, float>(fft_img, std::norm));
    Matrix<unsigned char> result(apply<std::complex<float>, float>(ifft_img, std::real));
    show_matrix(renderer, Matrix<unsigned char>(resultfft*255.0f/max(resultfft)));
    show_matrix(renderer, result);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(pWindow);
    SDL_Quit();
}

void test_filtrage_lineaire()
{
    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window* pWindow = SDL_CreateWindow("Image Processing Library",
                                            SDL_WINDOWPOS_UNDEFINED,
                                            SDL_WINDOWPOS_UNDEFINED,
                                            640,
                                            480,
                                            SDL_WINDOW_SHOWN);

    SDL_Renderer* renderer = SDL_CreateRenderer(pWindow, -1, SDL_RENDERER_ACCELERATED);

    std::size_t filtersize = 7;
    std::pair<Matrix<float>, Matrix<float> > XY = meshgrid<float>(filtersize);
    Matrix<float>& X = XY.first;
    Matrix<float>& Y = XY.second;
    X-=filtersize/2.0f-0.5f;
    Y-=filtersize/2.0f-0.5f;
    std::cout<<X<<Y<<std::endl;

    std::cout<<"savgol"<<std::endl;
    XY = meshgrid<float>(7,10);
    X = XY.first;
    Y = XY.second;
    Matrix<float> I = 5.0f+X*2.0f+Y*3.0f+7.0f*X*X+6.0f*Y*Y+9.0f*X*Y;
    std::cout<<I<<std::endl;
    std::cout<<xcorr(I,savgol(3,2,0,0))<<std::endl;
    std::cout<<2.0f+14.0f*X+9.0f*Y<<std::endl;
    std::cout<<xcorr(I,savgol(3,2,1,0))<<std::endl;
    std::cout<<3.0f+12.0f*Y+9.0f*X<<std::endl;
    std::cout<<xcorr(I,savgol(3,2,0,1))<<std::endl;
    std::cout<<14.0f*ones<float>(7,10)<<std::endl;
    std::cout<<xcorr(I,savgol(3,2,2,0))<<std::endl;
    std::cout<<9.0f*ones<float>(7,10)<<std::endl;
    std::cout<<xcorr(I,savgol(3,2,1,1))<<std::endl;
    std::cout<<12.0f*ones<float>(7,10)<<std::endl;
    std::cout<<xcorr(I,savgol(3,2,0,2))<<std::endl;

    // smoothing filters
    std::cout<<"smoothing filters"<<std::endl;
    std::cout<<savgol(7, 2)<<std::endl;
    std::cout<<average()<<std::endl;
    std::cout<<disk()<<std::endl;
    std::cout<<binomial()<<std::endl;
    std::cout<<gaussian()<<std::endl;
    std::cout<<pyramidal()<<std::endl;
    std::cout<<savgol(7,3)<<std::endl;

    std::cout<<"derivative filters x"<<std::endl;
    std::cout<<sobelx()<<std::endl;
    std::cout<<prewittx()<<std::endl;
    std::cout<<isotropicx()<<std::endl;
    std::cout<<mdifx()<<std::endl;
    std::cout<<savgol(7, 3, 1, 0)<<std::endl;

    std::cout<<"derivative filters y"<<std::endl;
    std::cout<<sobely()<<std::endl;
    std::cout<<prewitty()<<std::endl;
    std::cout<<isotropicy()<<std::endl;
    std::cout<<mdify()<<std::endl;
    std::cout<<savgol(7, 3, 0, 1)<<std::endl;

    std::cout<<"derivative filters xx"<<std::endl;
    std::cout<<mdifxx()<<std::endl;
    std::cout<<savgol(7, 3, 2, 0)<<std::endl;

    std::cout<<"derivative filters yy"<<std::endl;
    std::cout<<mdifyy()<<std::endl;
    std::cout<<savgol(7, 3, 0, 2)<<std::endl;

    std::cout<<"derivative filters xy"<<std::endl;
    std::cout<<mdifxy()<<std::endl;
    std::cout<<savgol(7, 3, 1, 1)<<std::endl;


    Matrix<Color> img = read_img("Images/Filtrage/phare_bruit_ps.png");
    Matrix<unsigned char> img_gray = color2grayimage(img);
    show_matrix(renderer, img_gray);
    show_matrix(renderer, filter(img_gray, average()));
    show_matrix(renderer, filter(img_gray, disk()));
    show_matrix(renderer, filter(img_gray, gaussian()));
    show_matrix(renderer, filter(img_gray, kirch()));
    show_matrix(renderer, filter(img_gray, laplacian(0.5)));
    show_matrix(renderer, filter(img_gray, log()));
    show_matrix(renderer, filter(img_gray, prewittx()));
    show_matrix(renderer, filter(img_gray, prewitty()));
    show_matrix(renderer, filter(img_gray, robinson()));
    show_matrix(renderer, filter(img_gray, sobelx()));
    show_matrix(renderer, filter(img_gray, sobely()));


    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(pWindow);
    SDL_Quit();
}

void test_filtrage_non_lineaire()
{
    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window* pWindow = SDL_CreateWindow("Image Processing Library",
                                            SDL_WINDOWPOS_UNDEFINED,
                                            SDL_WINDOWPOS_UNDEFINED,
                                            640,
                                            480,
                                            SDL_WINDOW_SHOWN);
    SDL_Renderer* renderer = SDL_CreateRenderer(pWindow, -1, SDL_RENDERER_ACCELERATED);

    Matrix<Color> img = read_img("Images/Filtrage/phare_bruit_ps.png");
    Matrix<unsigned char> img_gray = color2grayimage(img);
    show_matrix(renderer, img_gray);
    show_matrix(renderer, bilateral(img_gray, 10.0, 1.0));
    show_matrix(renderer, despeckle(img_gray));
    show_matrix(renderer, nagao(img_gray));

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(pWindow);
    SDL_Quit();
}

void test_snake()
{
    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window* pWindow = SDL_CreateWindow("Image Processing Library",
                                            SDL_WINDOWPOS_UNDEFINED,
                                            SDL_WINDOWPOS_UNDEFINED,
                                            640,
                                            480,
                                            SDL_WINDOW_SHOWN);

    SDL_Renderer* renderer = SDL_CreateRenderer(pWindow, -1, SDL_RENDERER_ACCELERATED);

    Matrix<Color> img = read_img("Images/Filtrage/phare_bruit_ps.png");
    Matrix<unsigned char> img_gray = color2grayimage(img);
    img_gray = nagao(img_gray);
    show_matrix(renderer, img_gray);
    Snake s;
    s.init_contour(150.0f, 150.0f, 50.0f);
    s.init_image(img_gray);
    std::size_t cpt=0;
    while(cpt++<20)
    {
        std::pair<Matrix<float>, Matrix<float> > Coordinates = s.getcoordinates();
        Matrix<float>& X = Coordinates.first, &Y = Coordinates.second;
        std::size_t n = Y.size();
        for(std::size_t i=0;i<n;++i)
            img(X(i, 0), Y(i, 0)) = Red;
        show_matrix(renderer, img);
        img = gray2colorimage(img_gray);
        s.iterate();
    }
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(pWindow);
    SDL_Quit();
}

void test_matrix()
{
    ///// tests construction de matrices
    Matrix<float> A = id<float>(3);
    std::cout<<A<<std::endl;
    A = zeros<float>(2, 4);
    std::cout<<A<<std::endl;
    A = ones<float>(7, 2);
    std::cout<<A<<std::endl;
    A = full<float>(4, 2.7f);
    std::cout<<A<<std::endl;
    A = arange<float>(1, 19);
    std::cout<<A<<std::endl;
    std::pair<Matrix<float>, Matrix<float> > XY = meshgrid<float>(3, 6);
    std::cout<<XY.first<<std::endl;
    std::cout<<XY.second<<std::endl;


    std::cout<<A<<std::endl;
    std::cout<<A.colNb()<<std::endl;
    std::cout<<A.rowNb()<<std::endl;
    std::cout<<A.size()<<std::endl;
    A.reshape(3,6);
    std::cout<<A<<std::endl;
    A.fliplr();
    std::cout<<A<<std::endl;
    A.flipud();
    std::cout<<A<<std::endl;
    A.rot90();
    std::cout<<A<<std::endl;
    A.rot180();
    std::cout<<A<<std::endl;
    A.rot270();
    std::cout<<A<<std::endl;
    A.transpose();
    std::cout<<A<<std::endl;
    A.setCol(1, full<float>(A.rowNb(), 2, 2.0));
    std::cout<<A<<std::endl;
}

void test_statistics()
{
    Matrix<float> M = arange(0.0f, 12.0f);
    M.reshape(3,4);
    std::cout<<M<<std::endl;
    std::cout<<min(M)<<std::endl;
    std::cout<<max(M)<<std::endl;
    std::cout<<mean(M)<<std::endl;
    std::cout<<var(M)<<std::endl;
    std::cout<<median(M)<<std::endl;
    std::cout<<geometric_mean(M)<<std::endl;
    std::cout<<sum(M)<<std::endl;
    std::cout<<prod(M)<<std::endl;
    std::cout<<cumsum(M)<<std::endl;
}

void test_linalg()
{
    Matrix<float> U = {{2,3,5},
                        {0,4,1},
                        {0,0,8}};

    Matrix<float> B = {{5,4,2,6,5},
                        {2,4,7,8,9},
                        {3,2,1,4,5}};

    Matrix<float> D;

    std::cout<<"pinv"<<std::endl;// pseudo inverse
    std::cout<<B<<std::endl;
    std::cout<<pinv(B)<<std::endl;
    std::cout<<dot(pinv(B), B)<<std::endl;
    std::cout<<dot(B, pinv(B))<<std::endl;

    std::cout<<"bwdsub"<<std::endl;// solve UX = B with U an upper triangular matrix
    Matrix<float> X = bwdsub(U, B);
    std::cout<<U<<std::endl;
    std::cout<<B<<std::endl;
    std::cout<<X<<std::endl;
    std::cout<<dot(U, X)<<std::endl;

    std::cout<<"fwdsub"<<std::endl;// solve LX = B with L a lower triangular matrix
    Matrix<float> L = transpose(U);
    X = fwdsub(L, B);
    std::cout<<L<<std::endl;
    std::cout<<B<<std::endl;
    std::cout<<X<<std::endl;
    std::cout<<dot(L, X)<<std::endl;

    std::cout<<"LU decomposition"<<std::endl;// LU decomposition PA = LU
    Matrix<float> A = {{2,3,5},
                        {6,4,1},
                        {2,9,8}};
    auto LUP = lu(A);
    Matrix<float> P;
    std::tie(L, U, P) = LUP;
    std::cout<<P<<std::endl;
    std::cout<<A<<std::endl;
    std::cout<<L<<std::endl;
    std::cout<<U<<std::endl;
    std::cout<<dot(P, A)<<std::endl;
    std::cout<<dot(L, U)<<std::endl;

    std::cout<<"QR decomposition"<<std::endl;// QR decomposition A = QR
    auto QR = qr(A);
    Matrix<float> Q,R;
    std::tie(Q,R) = QR;
    std::cout<<A<<std::endl;
    std::cout<<Q<<std::endl;
    std::cout<<dot(Q, transpose(Q))<<std::endl;
    std::cout<<dot(transpose(Q), Q)<<std::endl;
    std::cout<<R<<std::endl;
    std::cout<<dot(Q, R)<<std::endl;

    std::cout<<"RQ decomposition"<<std::endl;// RQ decomposition A = RQ
    auto RQ = rq(A);
    std::tie(R, Q) = RQ;
    std::cout<<A<<std::endl;
    std::cout<<Q<<std::endl;
    std::cout<<dot(Q, transpose(Q))<<std::endl;
    std::cout<<dot(transpose(Q), Q)<<std::endl;
    std::cout<<R<<std::endl;
    std::cout<<dot(R, Q)<<std::endl;

    std::cout<<"LQ decomposition"<<std::endl;// LQ decomposition A = LQ
    auto LQ = lq(A);
    std::tie(L, Q) = LQ;
    std::cout<<A<<std::endl;
    std::cout<<Q<<std::endl;
    std::cout<<dot(Q, transpose(Q))<<std::endl;
    std::cout<<dot(transpose(Q), Q)<<std::endl;
    std::cout<<L<<std::endl;
    std::cout<<dot(L, Q)<<std::endl;

    std::cout<<"QL decomposition"<<std::endl;// RQ decomposition A = QL
    auto QL = ql(A);
    std::tie(Q, L) = QL;
    std::cout<<A<<std::endl;
    std::cout<<Q<<std::endl;
    std::cout<<dot(Q, transpose(Q))<<std::endl;
    std::cout<<dot(transpose(Q), Q)<<std::endl;
    std::cout<<L<<std::endl;
    std::cout<<dot(Q, L)<<std::endl;

    std::cout<<"crout"<<std::endl;// Crout decomposition A = LDU
    auto LD = crout(dot(A, transpose(A)));
    std::tie(L, D) = LD;
    U = transpose(L);
    std::cout<<dot(A, transpose(A))<<std::endl;
    std::cout<<L<<std::endl;
    std::cout<<D<<std::endl;
    std::cout<<dot(dot(L, D),U)<<std::endl;


    std::cout<<"decomposition de cholesky"<<std::endl;// Cholesky decomposition A = LL'
    A = {{1, 2, 3},
        {2, 20, 26},
        {3, 26, 70}};
    L = cholesky(A);
    std::cout<<A<<std::endl;
    std::cout<<L<<std::endl;
    std::cout<<dot(L, transpose(L))<<std::endl;

    std::cout<<"jacobi"<<std::endl;// Jacobi decomposition A = LDL'
    A = {{3, 6, 10},
        {7, 3, 3},
        {7, 0, 0},
        {0, 3, 10},
        {9, 9, 2}};
    auto J = jacobi(dot(transpose(A),A));
    std::cout<<J.first<<std::endl;
    std::cout<<J.second<<std::endl;
    std::cout<<dot(transpose(A),A)<<std::endl;
    std::cout<<dot(J.second, dot(J.first, transpose(J.second)))<<std::endl;

    std::cout<<"bidiagonalisation"<<std::endl;
    auto ubv = bidiag_reduction(A);
    Matrix<float> V;
    std::tie(U,B,V) = ubv;
    std::cout<<U<<std::endl;
    std::cout<<B<<std::endl;
    std::cout<<V<<std::endl;
    std::cout<<A<<std::endl;
    std::cout<<dot(U, dot(B, V))<<std::endl;


    std::cout<<"svd"<<std::endl;
    auto USV = svd(A);
    Matrix<float> S;
    std::tie(U,S,V) = USV;
    std::cout<<A<<std::endl;
    std::cout<<U<<std::endl;
    std::cout<<S<<std::endl;
    std::cout<<V<<std::endl;
    std::cout<<dot(U, dot(S, transpose(V)))<<std::endl;

}

void test_regression()
{
    Matrix<float> X = arange<float>(-10.0,10.0);
    Matrix<float> Y = 10.0f*X*X+X*5.0f+2.0f;
    Y+=randn(X.rowNb(), X.colNb(), 0.0f, 0.3f);

    std::cout<<X<<std::endl;
    std::cout<<Y<<std::endl;
    std::cout<<poly_regression(Y, X, 2)<<std::endl;
}

void test_segmentation()
{
    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window* pWindow = SDL_CreateWindow("Image Processing Library",
                                            SDL_WINDOWPOS_UNDEFINED,
                                            SDL_WINDOWPOS_UNDEFINED,
                                            640,
                                            480,
                                            SDL_WINDOW_SHOWN);

    SDL_Renderer* renderer = SDL_CreateRenderer(pWindow, -1, SDL_RENDERER_ACCELERATED);


    Matrix<unsigned char> testimg = {{0,0,1,4,4,5},
                                    {0,1,3,4,3,4},
                                    {1,3,4,2,1,3},
                                    {4,4,3,1,0,0},
                                    {5,4,2,1,0,0},
                                    {5,5,4,3,1,0}};

    testimg*=51;
    testimg = read_pgm("Images/Pixmap/chat.pgm");
    testimg = color2grayimage(read_img("Images/Detection_de_forme/bruit1.bmp"));
    testimg = nagao(testimg);
    show_matrix(renderer, testimg);

    ///// test Otsu
    show_matrix(renderer, otsu(testimg));

    ///// test hysteresis
    show_matrix(renderer, hysteresis(testimg, 2*51, 3*51));

    // test felzenszwalb
    show_matrix(renderer, felzenszwalb(testimg));


    ///// test watershed
    Matrix<unsigned char> testwatershed = {{1,0,0,2,3,3,1,2,2,3,3}};
    Matrix<std::size_t>   markers       = {{0,1,1,0,0,0,2,0,0,0,0}};
    std::cout<<testwatershed<<std::endl;
    std::cout<<markers<<std::endl;
    std::cout<<watershed(testwatershed, markers)<<std::endl;


    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(pWindow);
    SDL_Quit();
}

void test_scaler()
{
    Matrix<unsigned char> testimg = {{0,0,1,4,4,5},
                                    {0,1,3,4,3,4},
                                    {1,3,4,2,1,3},
                                    {4,4,3,1,0,0},
                                    {5,4,2,1,0,0},
                                    {5,5,4,3,1,0}};

    // test scaler
    std::cout<<"test scaler"<<std::endl;
    MinMaxScaler sc(0, 255, 0);
    std::cout<<sc.fit_transform(Matrix<float>(testimg))<<std::endl;
    std::cout<<sc.inverse_transform(sc.fit_transform(Matrix<float>(testimg)))<<std::endl;
    sc = MinMaxScaler(0, 255, 1);
    std::cout<<sc.fit_transform(Matrix<float>(testimg))<<std::endl;
    std::cout<<sc.inverse_transform(sc.fit_transform(Matrix<float>(testimg)))<<std::endl;
    sc = MinMaxScaler(0, 255, 2);
    std::cout<<sc.fit_transform(Matrix<float>(testimg))<<std::endl;
    std::cout<<sc.inverse_transform(sc.fit_transform(Matrix<float>(testimg)))<<std::endl;
}

void test_hough()
{
    //test image reader
    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window* pWindow = SDL_CreateWindow("Image Processing Library",
                                            SDL_WINDOWPOS_UNDEFINED,
                                            SDL_WINDOWPOS_UNDEFINED,
                                            640,
                                            480,
                                            SDL_WINDOW_SHOWN);

    SDL_Renderer* renderer = SDL_CreateRenderer(pWindow, -1, SDL_RENDERER_ACCELERATED);




    //// test hough
    Matrix<bool> houghtest = {{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            };

    Matrix<bool> houghtest1 = {{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            };

    Matrix<bool> houghtest2 = {{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            };

    Matrix<bool> houghtest3 = {{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            };

    float phi_max = 2*(houghtest.rowNb()+houghtest.colNb())-1;

    Matrix<bool> normed = ones<bool>(houghtest.rowNb(),houghtest.colNb());
    Matrix<float> res2 = Matrix<float>(hough_transform(normed));

    Matrix<float> res = Matrix<float>(hough_transform(houghtest));
    std::transform(res.begin(), res.end(), res2.begin(), res.begin(), [](float v1, float v2){return v2>0.5?v1/v2:0.0f;});
    std::cout<<argmax(res)<<std::endl;
    std::cout<<Matrix<float>(argmax(res)).getCol(1)*2*PI/16.0f<<std::endl;
    std::cout<<Matrix<float>(argmax(res)).getCol(0)-phi_max/2.0f<<std::endl;
    show_matrix(renderer, houghtest);
    show_matrix(renderer, Matrix<unsigned char>(res*255/max(res)));
    show_matrix(renderer, Matrix<unsigned char>(res2));
    show_matrix(renderer, Matrix<unsigned char>(res*255.0/max(res)));
    show_matrix(renderer, label(Matrix<unsigned char>(res*255.0/max(res))));

    res = Matrix<float>(hough_transform(houghtest1));
    std::transform(res.begin(), res.end(), res2.begin(), res.begin(), [](float v1, float v2){return v2>0.5?v1/v2:0.0f;});
    std::cout<<argmax(res)<<std::endl;
    std::cout<<Matrix<float>(argmax(res)).getCol(1)*2*PI/16.0f<<std::endl;
    std::cout<<Matrix<float>(argmax(res)).getCol(0)-phi_max/2.0f<<std::endl;
    show_matrix(renderer, houghtest1);
    show_matrix(renderer, Matrix<unsigned char>(res*255.0/max(res)));

    res = Matrix<float>(hough_transform(houghtest2));
    std::transform(res.begin(), res.end(), res2.begin(), res.begin(), [](float v1, float v2){return v2>0.5?v1/v2:0.0f;});
    std::cout<<argmax(res)<<std::endl;
    std::cout<<Matrix<float>(argmax(res)).getCol(1)*2*PI/16.0f<<std::endl;
    std::cout<<Matrix<float>(argmax(res)).getCol(0)-phi_max/2.0f<<std::endl;
    show_matrix(renderer, houghtest2);
    show_matrix(renderer, Matrix<unsigned char>(res*255.0/max(res)));

    res = Matrix<float>(hough_transform(houghtest3));
    std::transform(res.begin(), res.end(), res2.begin(), res.begin(), [](float v1, float v2){return v2>0.5?v1/v2:0.0f;});
    std::cout<<argmax(res)<<std::endl;
    std::cout<<Matrix<float>(argmax(res)).getCol(1)*2*PI/16.0f<<std::endl;
    std::cout<<Matrix<float>(argmax(res)).getCol(0)-phi_max/2.0f<<std::endl;
    show_matrix(renderer, houghtest3);
    show_matrix(renderer, Matrix<unsigned char>(res*255.0/max(res)));

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(pWindow);
    SDL_Quit();
}

void test_zernike()
{
    //test image reader
    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window* pWindow = SDL_CreateWindow("Image Processing Library",
                                            SDL_WINDOWPOS_UNDEFINED,
                                            SDL_WINDOWPOS_UNDEFINED,
                                            640,
                                            480,
                                            SDL_WINDOW_SHOWN);

    SDL_Renderer* renderer = SDL_CreateRenderer(pWindow, -1, SDL_RENDERER_ACCELERATED);

    Matrix<bool> gray_img = {{0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,1,1,1,1,1,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0},};


    std::size_t H = gray_img.rowNb(), W = gray_img.colNb();

    show_matrix(renderer, gray_img);
    Matrix<std::size_t> indices = argwhere(gray_img);
    Matrix<float> gray_img_f(gray_img);
    Matrix<float> X(indices.getCol(0));
    Matrix<float> Y(indices.getCol(1));
    X-=float(H-1)/2.0f;
    Y-=float(W-1)/2.0f;
    X/=sqrt(float((H-1)*(H-1)+(W-1)*(W-1))/4.0f);
    Y/=sqrt(float((H-1)*(H-1)+(W-1)*(W-1))/4.0f);
    Matrix<float> rho = sqrt(X*X+Y*Y);
    Matrix<float> phi = atan2(Y, X);
    std::cout<<X<<" "<<Y<<std::endl;

    Matrix<std::complex<float> > F = zeros<std::complex<float> >(H, W);
    auto mesh = meshgrid<float>(H, W);
    Matrix<float>& meshX = mesh.first;
    Matrix<float>& meshY = mesh.second;
    meshX-=float(H-1)/2.0f;
    meshY-=float(W-1)/2.0f;
    meshX/=sqrt(float((H-1)*(H-1)+(W-1)*(W-1))/4.0f);
    meshY/=sqrt(float((H-1)*(H-1)+(W-1)*(W-1))/4.0f);
    Matrix<float> meshrho = sqrt(meshX*meshX+meshY*meshY);
    Matrix<float> meshphi = atan2(meshY, meshX);

    for(int i=0;i<6;++i)
    {
        for(int j=-i;j<=i;++j)
        {
            if((i-std::abs(j))%2==0)
            {
                Matrix<std::complex<float> > Aij(H, W);
                Zernike zer = Zernike(j, i);
                std::transform(rho.begin(), rho.end(), phi.begin(), Aij.begin(), [&zer](float r1, float r2){return zer.polynom(r1, r2);});
                Aij = apply<std::complex<float>, std::complex<float> >(Aij, std::conj);
                std::complex<float> aij = sum(Aij)*float(i+1)/float(PI);
                Matrix<std::complex<float> > f(H, W);
                std::transform(meshrho.begin(), meshrho.end(), meshphi.begin(), f.begin(), [&zer](float r1, float r2){return zer.polynom(r1, r2);});
                F+=f*aij;
            }
        }
        std::cout<<apply<std::complex<float>, float>(F, std::norm)<<std::endl;
    }

    Matrix<float> F_float = apply<std::complex<float>, float>(F, std::norm);
    F_float/=max(F_float)/255.0f;
    show_matrix(renderer, Matrix<unsigned char>(F_float));
    std::cout<<apply<std::complex<float>, float>(F, std::norm)<<std::endl;
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(pWindow);
    SDL_Quit();
}

void test_digits()
{
    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window* pWindow = SDL_CreateWindow("Image Processing Library",
                                            SDL_WINDOWPOS_UNDEFINED,
                                            SDL_WINDOWPOS_UNDEFINED,
                                            480,
                                            480,
                                            SDL_WINDOW_SHOWN);

    SDL_Renderer* renderer = SDL_CreateRenderer(pWindow, -1, SDL_RENDERER_ACCELERATED);

   std::vector<std::string> filenames = get_files_recursively("Images/Digits/", ".png");
   std::size_t nb_samples = filenames.size();
   Matrix<float> Data(nb_samples, 64);
   Matrix<std::size_t> Labels(nb_samples, 1);
   for(std::size_t i=0;i<nb_samples;++i)
   {
       std::cout<<filenames[i]<<" "<<filenames[i][14]<<std::endl;
       Matrix<Color> img = read_png(filenames[i]);
       //show_matrix(renderer, img);
       Matrix<float> gray_img(color2grayimage(img));
       gray_img.reshape(1, 64);
       Data.setRow(i, gray_img/127.5f-1.0f);
       Labels(i, 0) = filenames[i][14]-'0';
   }
    MLP clf({});
    Matrix<std::size_t> Result = clf.fit_predict(Data, Labels);
    clf.export_graphviz("digits.dot");
    std::cout<<count(Result==Labels, true)<<"/"<<nb_samples<<std::endl;
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(pWindow);
    SDL_Quit();
}

void test_polynomials()
{
    std::cout<<Polynomial({0.0f, 2.0f})<<std::endl;

    for(unsigned i=0;i<10;++i)
        std::cout<<Tchebychev(i)<<std::endl;
    for(unsigned i=0;i<10;++i)
        std::cout<<Laguerre(i)<<std::endl;
    for(unsigned i=0;i<10;++i)
        std::cout<<Legendre(i)<<std::endl;
    for(unsigned i=0;i<10;++i)
        std::cout<<Hermite_pro(i)<<std::endl;
    for(unsigned i=0;i<10;++i)
    {
        std::cout<<Hermite_phi(i).integrate()<<std::endl;
        std::cout<<Hermite_phi(i)<<std::endl;
        std::cout<<Hermite_phi(i).differentiate()<<std::endl;
    }

}

void test_pca()
{
    PCA pca(1);
    Matrix<float> X = arange<float>(-10.0,10.0);
    Matrix<float> Y = 5.0f*X+2.0f;
    Y+=randn(Y.rowNb(), Y.colNb(), 0.0f, 0.3f);
    Matrix<float> Data = horzcat(X,Y);
    Matrix<float> Z = pca.fit_transform(Data);

    std::cout<<Data<<std::endl;
    std::cout<<Z<<std::endl;
    std::cout<<pca.inverse_transform(Z)<<std::endl;
}

void test_nnmf()
{
    NNMF nnmf(1);
    Matrix<float> X = arange<float>(0.0,100.0);
    Matrix<float> Y = 5.0f*X+2.0f;
    Y+=randn(Y.rowNb(), Y.colNb(), 0.0f, 0.3f);
    Matrix<float> Data = horzcat(X,Y);
    Matrix<float> Z = nnmf.fit_transform(Data);

    std::cout<<Data<<std::endl;
    std::cout<<Z<<std::endl;
    std::cout<<nnmf.inverse_transform(Z)<<std::endl;
}

void test_linprog()
{
    Matrix<float> f = {{-1},{-2}};
    Matrix<float> A = {{-3,2},{-1,2},{1,1}};
    Matrix<float> b = {{2},{4},{5}};
    std::cout<<linprog_can(f, A, b)<<std::endl;// answer {{2},{3}}

    f = {{-3},{-2}};
    A = {{2,1},{1,1},{1,0}};
    b = {{100},{80},{40}};
    std::cout<<linprog_can(f, A, b)<<std::endl;// answer {{20},{60}}

    f = {{-6},{-4}};
    A = {{3,9},{4,5},{2,1}};
    b = {{81},{55},{20}};
    std::cout<<linprog_can(f, A, b)<<std::endl;// answer {{15/2},{5}}

    f = {{1},{1}};
    A = {{-1,-2}};
    b = {{-2}};
    std::cout<<linprog_can(f, A, b)<<std::endl;// answer {{0},{1}}

    f = {{2},{-5}};
    A = {{1,0},{-1,0},{0,1},{0,-1},{-1,-1}};
    b = {{200},{-100},{170},{-80},{-200}};
    std::cout<<linprog_can(f, A, b)<<std::endl;// answer {{100},{170}}
}

void test_quadprog()
{
    Matrix<float> H = {{2,0},{0,2}};
    Matrix<float> f = {{-2},{-4}};
    Matrix<float> A = {{1,2}};
    Matrix<float> b = {{4}};
    Matrix<float> Aeq = {{0,0}};
    Matrix<float> beq = {{0}};
    Matrix<float> x0 = {{2},{0}};
    std::cout<<-dot(pinv(H), f)<<std::endl;
    std::cout<<conjugate_gradient(H,f,x0)<<std::endl;
    std::cout<<solve(H,f,x0)<<std::endl;
    std::cout<<quadprog(H, f, A, b, Aeq, beq, x0)<<std::endl;

    H = {{3.0f,1.0f},{1.0f,8.0f}};
    f = {{42.6f},{2.7f}};
    x0 = {{34.2f},{34.2f}};
    std::cout<<conjugate_gradient(H,f,x0)<<std::endl;//14.7 -1.5
    std::cout<<solve(H,f,x0)<<std::endl;//14.7 -1.5
    std::cout<<dot(A,conjugate_gradient(H,f,x0))<<std::endl;
}

void test_interp()
{
    Matrix<float> X = {{0.0f,1.0f,2.0f,3.0f,4.0f,5.0f,6.0f,7.0f,8.0f,9.0f,10.0f}};
    Matrix<float> Y = X*4.2f+1.3f;
    Matrix<float> Z = 2.0f*X*X+5.0f*X+3.0f;
    Matrix<float> X2 = X+0.4f;
    X2.newCol(0);
    X2(0, 0) = -0.6f;
    Matrix<float> Y2 = X2*4.2f+1.3f;
    Matrix<float> Z2 = 2.0f*X2*X2+5.0f*X2+3.0f;
    std::cout<<X<<std::endl;
    std::cout<<X2<<std::endl;
    std::cout<<Y<<std::endl;
    std::cout<<Y2<<std::endl;
    std::cout<<interp1(X, Y, X2, "linear")<<std::endl;
    std::cout<<interp1(X, Y, X2, "nearest")<<std::endl;
    std::cout<<interp1(X, Y, X2, "cubic")<<std::endl;
    std::cout<<Z<<std::endl;
    std::cout<<Z2<<std::endl;
    std::cout<<interp1(X, Z, X2, "linear")<<std::endl;
    std::cout<<interp1(X, Z, X2, "nearest")<<std::endl;
    std::cout<<interp1(X, Z, X2, "cubic")<<std::endl;

}


void test_optim()
{
    Matrix<float> X0 = {{34.2f},{34.2f}};
    std::cout<<"himmelblau"<<std::endl;
    std::cout<<steepest_descent(himmelblau, X0)<<std::endl;
    std::cout<<"rosenbrock"<<std::endl;
    std::cout<<steepest_descent(rosenbrock, X0)<<std::endl;
    std::cout<<"himmelblau lbfgs"<<std::endl;
    std::cout<<lbfgs(himmelblau, X0)<<std::endl;
    std::cout<<"rosenbrock lbfgs"<<std::endl;
    std::cout<<lbfgs(rosenbrock, X0)<<std::endl;
    std::cout<<"himmelblau bfgs"<<std::endl;
    std::cout<<bfgs(himmelblau, X0)<<std::endl;
    std::cout<<"rosenbrock bfgs"<<std::endl;
    std::cout<<bfgs(rosenbrock, X0)<<std::endl;
}
