#define SDL_MAIN_HANDLED
#include <SDL2/SDL.h>
#include <cfloat>

#include "csvIO.hpp"
#include "file_explorer.hpp"
#include "linalg.hpp"
#include "logical_operators.hpp"
#include "matrix.hpp"
#include "mmath.hpp"
#include "polynomial_equation_solver.hpp"
#include "regression.hpp"
#include "statistics.hpp"

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
#include "ImageProcessing/imgconverter.hpp"
#include "ImageProcessing/linear_filtering.hpp"
#include "ImageProcessing/nonlinear_filtering.hpp"
#include "ImageProcessing/morphology.hpp"
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


char fctbidon(float a)
{
    char c = a>3?'a':'d';
    return c;
}

char fctbidon(bool a)
{
    char c = a?'a':'d';
    return c;
}


SDL_Window* init_SDL()
{
    SDL_Init(SDL_INIT_VIDEO);
    return SDL_CreateWindow("hello", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, 800, 600, SDL_WINDOW_SHOWN);
}

SDL_Renderer* get_SDL_renderer(SDL_Window* screen)
{
    return SDL_CreateRenderer(screen, -1, SDL_RENDERER_PRESENTVSYNC);
}

void pause()
{
    bool done = false;
    SDL_Event event;
    while(!done)
    {
        while(SDL_PollEvent(&event))
        {
            if(event.type==SDL_QUIT)
                done = true;
        }
    }
}

void show_matrix(SDL_Renderer* renderer, const Matrix<Color>& img)
{
    SDL_RenderClear(renderer);
    SDL_Texture* tex = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ABGR8888, SDL_TEXTUREACCESS_STREAMING, img.colNb(), img.rowNb());
    Color* pColor = NULL;

    pColor = (Color*) malloc(sizeof(Color)*img.size());
    std::copy(img.cbegin(), img.cend(), pColor);

    SDL_UpdateTexture(tex, NULL, pColor, img.colNb()*sizeof(Color));
    SDL_RenderCopy(renderer, tex, NULL, NULL);
    SDL_RenderPresent(renderer);
    SDL_DestroyTexture(tex);
    free(pColor);
    pause();
}



void show_matrix(SDL_Renderer* renderer, const Matrix<unsigned char>& img)
{
    show_matrix(renderer, gray2colorimage(img));
}

void show_matrix(SDL_Renderer* renderer, const Matrix<bool>& img)
{
    show_matrix(renderer, bw2colorimage(img));
}

void show_matrix(SDL_Renderer* renderer, const Matrix<int>& img)
{
    show_matrix(renderer, array2colorimage(img));
}


int main()
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

    //// test filter
    Matrix<unsigned char> A_img = rand<unsigned char>(3,4);
    std::cout<<Matrix<int>(A_img)<<std::endl;
    Matrix<float> filter = zeros<float>(3,3);
    filter(0,0)=1.0f;
    Matrix<float> Af_img = conv(Matrix<float>(A_img), filter);
    std::cout<<Af_img<<std::endl;

    ///////////test clustering

    //create data
    int n_samples = 100;
    int n_features = 2;
    int n_classes = 4;
    Matrix<float> Data(n_samples, n_features);
    Matrix<int> Labels(n_samples, 1);
    Matrix<float> means = {{10.0f, 20.0f, 30.0f, 40.0f},
                            {15.0f, 40.0f, 5.0f, 30.0f}};
    Matrix<float> stddevs = {{1.0f, 2.0f, 3.0f, 4.0f},
                            {2.5f, 1.5f, 1.0f, 2.0f}};

    Matrix<float> results;

    /*Matrix<float> means = {{10.0f, 20.0f},
                            {15.0f, 40.0f}};
    Matrix<float> stddevs = {{1.0f, 2.0f},
                            {2.5f, 1.5f}};*/

    stddevs = full<float>(n_features, n_classes, 1.0f);
    for(int i=0;i<n_classes;++i)
    {
        for(int j=0;j<n_features;++j)
            Data.setSubmat(i*n_samples/n_classes, j, randn(n_samples/n_classes, 1, means(j, i), stddevs(j, i)));
        Labels.setSubmat(i*n_samples/n_classes, 0, full<int>(n_samples/n_classes, 1, i));
    }


    // test Affinity Propagation
    std::cout<<"Affinity Propagation"<<std::endl;
    AffinityPropagation clf_affpro;
    results = Matrix<float>(clf_affpro.fit_predict(Data));
    results.reshape(n_classes, n_samples/n_classes);
    std::cout<<histogram(Matrix<unsigned char>(results)).getCols(0, n_classes)<<std::endl;
    std::cout<<axismean(results, 1)<<std::endl;
    std::cout<<axisstdev(results, 1)<<std::endl;
    std::cout<<clf_affpro.getCenters()<<std::endl;

    // test Kmeans
    std::cout<<"Kmeans"<<std::endl;
    Kmeans clf(n_classes);
    results = Matrix<float>(clf.fit_predict(Data));
    results.reshape(n_classes, n_samples/n_classes);
    std::cout<<histogram(Matrix<unsigned char>(results)).getCols(0, n_classes)<<std::endl;
    std::cout<<axismean(results, 1)<<std::endl;
    std::cout<<axisstdev(results, 1)<<std::endl;
    std::cout<<clf.getCenters()<<std::endl;

    // test, Kmedoids
    std::cout<<"Kmedoids"<<std::endl;
    Kmedoids kmed_clf(n_classes);
    results = Matrix<float>(kmed_clf.fit_predict(Data));
    results.reshape(n_classes, n_samples/n_classes);
    std::cout<<histogram(Matrix<unsigned char>(results)).getCols(0, n_classes)<<std::endl;
    std::cout<<axismean(results, 1)<<std::endl;
    std::cout<<axisstdev(results, 1)<<std::endl;
    std::cout<<kmed_clf.getCenters()<<std::endl;

    // test hierarchical clustering
    std::cout<<"hierarchical clustering"<<std::endl;
    HierarchicalClustering hier_clf(n_classes);
    results = Matrix<float>(hier_clf.fit_predict((Data)));
    results.reshape(n_classes, n_samples/n_classes);
    std::cout<<histogram(Matrix<unsigned char>(results)).getCols(0, n_classes)<<std::endl;
    std::cout<<axismean(results, 1)<<std::endl;
    std::cout<<axisstdev(results, 1)<<std::endl;

    //test gmm
    std::cout<<"GMM"<<std::endl;
    GMM clf_gmm(n_classes);
    results = Matrix<float>(clf_gmm.fit_predict(Data));
    results.reshape(n_classes, n_samples/n_classes);
    std::cout<<histogram(Matrix<unsigned char>(results)).getCols(0, n_classes)<<std::endl;
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
    std::cout<<histogram(Matrix<unsigned char>(results)).getCols(0, n_classes)<<std::endl;
    std::cout<<axismean(results, 1)<<std::endl;
    std::cout<<axisstdev(results, 1)<<std::endl;



    ///// test gaussian_naives_bayes
    std::cout<<"Gaussian naives bayes"<<std::endl;
    GaussianNaiveBayes clf_gnb;
    clf_gnb.fit(Data, Labels);
    results = Matrix<float>(clf_gnb.fit_predict(Data, Labels));
    std::cout<<count_nonzero(Matrix<int>(results) == Labels)<<std::endl;
    std::cout<<clf_gnb.getCenters()<<std::endl;
    std::cout<<clf_gnb.getStddevs()<<std::endl;
    std::cout<<clf_gnb.getPriors()<<std::endl;

    ///// test Dtree
    std::cout<<"Dtree"<<std::endl;
    Dtree tree;
    results = Matrix<float>(tree.fit_predict(Data, Labels));
    std::cout<<count_nonzero(Matrix<int>(results) == Labels)<<std::endl;
    tree.export_graphviz("tree.dot");

    ///// test k-NN
    std::cout<<"KNN"<<std::endl;
    KNN clf_knn(5);
    clf_knn.fit(Data, Labels);
    results = Matrix<float>(clf_knn.fit_predict(Data, Labels));
    std::cout<<count_nonzero(Matrix<int>(results) == Labels)<<std::endl;

    ///// test MLP
    std::cout<<"MLP"<<std::endl;
    StandardScaler ssc(1,1, 2);
    Data = ssc.fit_transform(Data);
    std::vector<int> param = {};
    MLP clf_mlp(param);
    results = Matrix<float>(clf_mlp.fit_predict(Data, Labels));
    clf_mlp.export_graphviz("test.dot");
    std::cout<<count_nonzero(Matrix<int>(results) == Labels)<<std::endl;

    return 0;

    ///// test Otsu
    Matrix<unsigned char> testimg = {{0,0,1,4,4,5},
                                    {0,1,3,4,3,4},
                                    {1,3,4,2,1,3},
                                    {4,4,3,1,0,0},
                                    {5,4,2,1,0,0},
                                    {5,5,4,3,1,0}};
    Matrix<Color> testimg2 = gray2colorimage(testimg);
    std::cout<<testimg2<<std::endl;
    std::cout<<otsu(testimg)<<std::endl;
    std::cout<< (testimg>(unsigned char)(4)) <<std::endl;
    std::cout<<argwhere(testimg>4)<<std::endl;
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

    //test image reader
    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window* pWindow = SDL_CreateWindow("Image Processing Library",
                                            SDL_WINDOWPOS_UNDEFINED,
                                            SDL_WINDOWPOS_UNDEFINED,
                                            640,
                                            480,
                                            SDL_WINDOW_SHOWN);

    SDL_Renderer* renderer = SDL_CreateRenderer(pWindow, -1, SDL_RENDERER_ACCELERATED);

    std::string extensions[] = {"bmp"};//, "bmp", "tga", "ico", "png", "pbm,"pgm", "ppm", "jpg", "gif", "tiff"};

    std::for_each(
                  extensions,
                  extensions+1,
                  [renderer](std::string s){auto cur_list = get_files_recursively(".", s);
                                    std::for_each(cur_list.rbegin()+12,
                                                  cur_list.rend(),
                                                  [renderer](std::string ss){std::cout<<ss<<std::endl;
                                                                            Matrix<Color> img = read_bmp(ss);
                                                                            //img = img.getSubmat(200, 300, 200, 300);
                                                                            show_matrix(renderer, img);
                                                                            Matrix<float> imgR = Matrix<float>(apply(img, &Color::red));
                                                                            Matrix<float> imgG = Matrix<float>(apply(img, &Color::green));
                                                                            Matrix<float> imgB = Matrix<float>(apply(img, &Color::blue));
                                                                            int n = img.size();
                                                                            Matrix<float> img_data(n, 3);
                                                                            imgR.reshape(n, 1);imgG.reshape(n, 1);imgB.reshape(n, 1);
                                                                            img_data.setCol(0, imgR);
                                                                            img_data.setCol(1, imgG);
                                                                            img_data.setCol(2, imgB);
                                                                            std::cout<<"start"<<std::endl;
                                                                            Kmeans kmeans_clf(3);
                                                                            Matrix<int> labels = kmeans_clf.fit_predict(img_data);
                                                                            std::cout<<histogram(labels)<<std::endl;
                                                                            std::cout<<"finished"<<std::endl;
                                                                            labels.reshape(img.rowNb(), img.colNb());
                                                                            Matrix<float> centers = kmeans_clf.getCenters();
                                                                            Matrix<unsigned char> centers_u = Matrix<unsigned char>(Matrix<int>(centers));
                                                                            std::cout<<centers_u<<std::endl;
                                                                            Matrix<Color> new_img(img.size(), 1);
                                                                            for(int a=0;a<img.rowNb();++a)
                                                                            {
                                                                                    int lab = labels(a, 0);
                                                                                    new_img(a, 0) = Color(centers_u(lab, 0), centers_u(lab, 1), centers_u(lab, 2));
                                                                            }
                                                                            new_img.reshape(img.rowNb(), img.colNb());
                                                                            labels.reshape(img.rowNb(), img.colNb());
                                                                            show_matrix(renderer, labels);
                                                                            show_matrix(renderer, new_img);
                                                                            }
                                                  );}
                  );





    //Test writer
    Matrix<Color> img = read_png("Images/Detection_de_forme/bild4.png");
    show_matrix(renderer, img);
    save_ppm("testwriter.ppm", img);
    save_pgm("testwriter.pgm", color2grayimage(img));
    save_pbm("testwriter.pbm", color2bwimage(img, 127));
    save_png("testwriter.png", img);
    save_jpeg("testwriter.jpg", img);
    save_tga("testwriter.tga", img);
    save_bmp("testwriter.bmp", img);

    // test filtres non linéaires
    img = read_png("Images/Filtrage/phare_bruit_ps.png");
    Matrix<unsigned char> gray_img = color2grayimage(img);
    show_matrix(renderer, gray_img);
    show_matrix(renderer, bilateral(gray_img, 5, 1000, 5));
    show_matrix(renderer, despeckle(gray_img));
    show_matrix(renderer, nagao(gray_img));
    show_matrix(renderer, conservative_smoothing(gray_img));
    //gray_img = opening_by_reconstruction(gray_img);
    //gray_img = closing_by_reconstruction(gray_img);

    // test binarisation
    Matrix<unsigned char> G = gradient(gray_img);
    show_matrix(renderer, G);
    Matrix<bool> Gb = hysteresis(G, 30, 100);
    show_matrix(renderer, Gb);
    show_matrix(renderer, otsu(gray_img));

    // test FFT
    Matrix<float> data = sin(2.0f*float(PI)/10.0f*arange<float>(0.0f,  10.0f));
    std::cout<<data<<std::endl;
    std::cout<<FFT(data)<<std::endl;
    std::cout<<iFFT(FFT(data))<<std::endl;


    //essai webcam
    /*SDL_Surface* image = SDL_CreateRGBSurface(SDL_SWSURFACE, 640, 480, 32, 0xFF000000, 0x00FF0000, 0x0000FF00, 0x000000FF);
    SDL_SysWMinfo wmInfo;
    SDL_VERSION(&wmInfo.version);
    SDL_GetWindowWMInfo(pWindow, &wmInfo);
    if(wmInfo.subsystem==SDL_SYSWM_WINDOWS)
        std::cout<<"sub ok"<<std::endl;
    auto hWnd = wmInfo.info.win.window;
    auto hDC = GetDC(hWnd);
    auto hWnd_WC = capCreateCaptureWindow("WebCam", WS_CHILD, 0, 0, 1, 1, hWnd, 0);

    if(!capDriverConnect(hWnd_WC, 0))
    {
        MessageBox(NULL, "Erreur", "Erreur", MB_ICONERROR);
        return FALSE;
    }*/



    SDL_DestroyWindow(pWindow);
    SDL_Quit();
    return 0;
}
