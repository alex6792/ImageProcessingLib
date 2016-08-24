#include "csvIO.hpp"
#include "file_explorer.hpp"
#include "linalg.hpp"
#include "logical_operators.hpp"
#include "matrix.hpp"
#include "mmath.hpp"
#include "polynomial_equation_solver.hpp"
#include "regression.hpp"
#include "SDL_interface.hpp"
#include "statistics.hpp"
#include "test_functions.hpp"

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
#include "MachineLearning/SupervisedLearning/mlp.hpp"
#include "MachineLearning/SupervisedLearning/perceptron.hpp"

#include "ImageProcessing/color.hpp"
#include "ImageProcessing/fourier.hpp"
#include "ImageProcessing/hough.hpp"
#include "ImageProcessing/imgconverter.hpp"
#include "ImageProcessing/linear_filtering.hpp"
#include "ImageProcessing/morphology.hpp"
#include "ImageProcessing/nonlinear_filtering.hpp"
#include "ImageProcessing/segmentation.hpp"
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


void test_clustering()
{
     ///////////test clustering
    //create data
    std::size_t n_samples = 400;
    std::size_t n_features = 2;
    std::size_t n_classes = 4;
    Matrix<float> Data(n_samples, n_features);
    Matrix<std::size_t> Labels(n_samples, 1);
    Matrix<float> means = {{10.0f, 20.0f, 30.0f, 40.0f},
                            {15.0f, 40.0f, 5.0f, 30.0f}};
    Matrix<float> stddevs = {{1.0f, 2.0f, 3.0f, 4.0f},
                            {2.5f, 1.5f, 1.0f, 2.0f}};

    Matrix<float> results;

    stddevs = full<float>(n_features, n_classes, 1.0f);
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

    //test gmm
    std::cout<<"GMM"<<std::endl;
    GMM clf_gmm(n_classes);
    results = Matrix<float>(clf_gmm.fit_predict(Data));
    results.reshape(n_classes, n_samples/n_classes);
        std::cout<<(histogram(Matrix<std::size_t>(results))).getCols(0, n_classes)<<std::endl;
    std::cout<<axismean(results, 1)<<std::endl;
    std::cout<<axisstdev(results, 1)<<std::endl;
    std::cout<<clf_gmm.getCenters()<<std::endl;
    std::cout<<clf_gmm.getStddevs()<<std::endl;
    std::cout<<clf_gmm.getWeights()<<std::endl;

    //test meanshift
    std::cout<<"Mean Shift"<<std::endl;
    MeanShift clf_meanshift = MeanShift();
    results = Matrix<float>(clf_meanshift.fit_predict(Data));
    results.reshape(n_classes, n_samples/n_classes);
        std::cout<<(histogram(Matrix<std::size_t>(results))).getCols(0, n_classes)<<std::endl;
    std::cout<<axismean(results, 1)<<std::endl;
    std::cout<<axisstdev(results, 1)<<std::endl;
}

void test_supervisedlearning()
{
        //create data
    std::size_t n_samples = 200;
    std::size_t n_features = 2;
    std::size_t n_classes = 4;
    Matrix<float> Data(n_samples, n_features);
    Matrix<std::size_t> Labels(n_samples, 1);
    Matrix<float> means = {{10.0f, 20.0f, 30.0f, 40.0f},
                            {15.0f, 40.0f, 5.0f, 30.0f}};
    Matrix<float> stddevs = {{1.0f, 2.0f, 3.0f, 4.0f},
                            {2.5f, 1.5f, 1.0f, 2.0f}};

    Matrix<float> results;

    stddevs = full<float>(n_features, n_classes, 1.0f);
    for(std::size_t i=0;i<n_classes;++i)
    {
        for(std::size_t j=0;j<n_features;++j)
            Data.setSubmat(i*n_samples/n_classes, j, randn(n_samples/n_classes, 1, means(j, i), stddevs(j, i)));
        Labels.setSubmat(i*n_samples/n_classes, 0, full<std::size_t>(n_samples/n_classes, 1, i));
    }

    ///// test gaussian_naives_bayes
    std::cout<<"Gaussian naives bayes"<<std::endl;
    GaussianNaiveBayes clf_gnb;
    clf_gnb.fit(Data, Labels);
    results = Matrix<float>(clf_gnb.fit_predict(Data, Labels));
    std::cout<<count_nonzero(Matrix<std::size_t>(results) == Labels)<<std::endl;
    std::cout<<clf_gnb.getCenters()<<std::endl;
    std::cout<<clf_gnb.getStddevs()<<std::endl;
    std::cout<<clf_gnb.getPriors()<<std::endl;

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
}

void test_morpho_binaire()
{
    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window* pWindow = SDL_CreateWindow("Image Processing Library",
                                            SDL_WINDOWPOS_UNDEFINED,
                                            SDL_WINDOWPOS_UNDEFINED,
                                            640,
                                            480,
                                            SDL_WINDOW_SHOWN);

    SDL_Renderer* renderer = SDL_CreateRenderer(pWindow, -1, SDL_RENDERER_ACCELERATED);

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
    Matrix<Color> img = read_png("Images/Filtrage/phare_bruit_ps.png");
    Matrix<unsigned char> gray_img = color2grayimage(img);

    SDL_SetWindowTitle(pWindow, "original");
    show_matrix(renderer, gray_img);
    SDL_SetWindowTitle(pWindow, "bilateral");
    show_matrix(renderer, bilateral(gray_img, 5, 1000, 5));
    SDL_SetWindowTitle(pWindow, "despeckle");
    show_matrix(renderer, despeckle(gray_img));
    SDL_SetWindowTitle(pWindow, "nagao");
    show_matrix(renderer, nagao(gray_img));
    SDL_SetWindowTitle(pWindow, "erode");
    show_matrix(renderer, erode(gray_img));
    SDL_SetWindowTitle(pWindow, "dilate");
    show_matrix(renderer, dilate(gray_img));
    SDL_SetWindowTitle(pWindow, "median");
    show_matrix(renderer, median_filter(gray_img));
    SDL_SetWindowTitle(pWindow, "conservative smoothing");
    show_matrix(renderer, conservative_smoothing(gray_img));
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

    std::string extensions[] = { "bmp", "tga", "ico", "png", "pbm","pgm", "ppm", "jpg", "gif", "tiff"};

    std::for_each(
                  extensions,
                  extensions+1,
                  [renderer](std::string s){auto cur_list = get_files_recursively(".", s);
                                                std::for_each(cur_list.begin(),
                                                  cur_list.end(),
                                                  [renderer](std::string ss){std::cout<<ss<<std::endl;
                                                                            Matrix<Color> img = read_img(ss);
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
    //Test writer
    Matrix<Color> img = read_png("Images/Detection_de_forme/bild4.png");
    std::vector<std::string> extensions = { "bmp", "tga", "ico", "png", "pbm","pgm", "ppm", "jpg", "gif", "tiff"};
    std::for_each(extensions.cbegin(), extensions.cend(), [img](const std::string& s){save_img("testwriter."+s, img);});

}

void test_fourier()
{
    // test FFT
    Matrix<float> data = sin(2.0f*float(PI)/10.0f*arange<float>(0.0f,  10.0f));
    std::cout<<data<<std::endl;
    std::cout<<FFT(data)<<std::endl;
    std::cout<<iFFT(FFT(data))<<std::endl;
}

void test_filtrage__lineaire()
{
    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window* pWindow = SDL_CreateWindow("Image Processing Library",
                                            SDL_WINDOWPOS_UNDEFINED,
                                            SDL_WINDOWPOS_UNDEFINED,
                                            640,
                                            480,
                                            SDL_WINDOW_SHOWN);

    SDL_Renderer* renderer = SDL_CreateRenderer(pWindow, -1, SDL_RENDERER_ACCELERATED);



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

}

void test_linalg()
{

}

void test_regression()
{

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

void test_sclaler()
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
}
